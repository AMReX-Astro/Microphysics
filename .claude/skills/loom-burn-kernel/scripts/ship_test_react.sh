#!/bin/bash
# Ship the local test_react Loom integration to the GPU host's Microphysics
# tree, show what differs, and optionally apply it and start the HIP build.
#
#   ship_test_react.sh [--apply] [--build] KERNEL_SUBDIR
#       Stage the local test_react files (GNUmakefile Make.package _parameters
#       main.cpp loom_*.H loom_validate.cpp) and KERNEL_SUBDIR (relative to
#       test_react, e.g. loom/<port>; its tests/ and build dirs are left out)
#       into $REMOTE_SHIP_ROOT/ship-<timestamp>-<pid> on the host and list the files
#       that are NEW or DIFFER from $REMOTE_TREE/unit_test/test_react.
#       --apply  also copy the NEW/DIFFER files into the remote tree.
#       --build  (implies --apply) then start the HIP build of test_react,
#                detached, logged to <ship dir>/build.log:
#                  make -j$JOBS USE_HIP=TRUE AMD_ARCH=$ARCH NETWORK_DIR=...
#                    EOS_DIR=... INTEGRATOR_DIR=... USE_NEUTRINOS=FALSE
#                    XTRA_CXXFLAGS=-ffp-contract=off USE_LOOM_BURN=TRUE
#                    LOOM_KERNEL_DIR=KERNEL_SUBDIR $EXTRA_MAKE
#   ship_test_react.sh status
#       Show whether the last --build is still running, and its log tail.
#   ship_test_react.sh -h
#
# With LOOM_HOST=local (the default) the "remote" tree is this checkout:
# nothing differs, and --build runs the HIP build in place, detached, with the
# same log and status handling.
#
# Remote: every invocation makes exactly ONE ssh connection: the payload goes through
# ssh's stdin and is unpacked in the foreground; only the build job is
# backgrounded (nohup setsid, stdin from /dev/null). Space status checks out.
#
# Environment:
#   LOOM_HOST         ssh alias of the GPU host, or "local"      (default local)
#   LOOM_SSH_OPTS     extra ssh options, word-split  (default: -o WarnWeakCrypto=no
#                     if the local OpenSSH is >= 10.1, else none)
#   LOOM_REMOTE_BASE  base dir on the GPU host (absolute)        (default $LOOM_WORK or $HOME/loom-work when
#                     local; required for a remote host, e.g. /mnt/data/bwibking on moth)
#   REMOTE_TREE       Microphysics checkout on the GPU host      (default: local: the git toplevel
#                                                                 of cwd; remote: $LOOM_REMOTE_BASE/Microphysics)
#   REMOTE_SHIP_ROOT  where ship-<timestamp>/ dirs are created   (default $LOOM_REMOTE_BASE/ship)
#   LOCAL_TEST_DIR    local test_react dir   (default <git toplevel of cwd>/unit_test/test_react)
#   --build only:
#   NETWORK_DIR, EOS_DIR, INTEGRATOR_DIR   make variables         (required, no defaults)
#   ARCH              AMD_ARCH        (default: the host's GPU from amdgpu-arch/rocminfo, else gfx90a)
#   JOBS              make -j                                    (default: nproc on the host)
#   ROCM_PATH         ROCm on the GPU host                       (default /opt/rocm)
#   AMREX_HOME        AMReX checkout on the GPU host             (default: local: $AMREX_HOME;
#                                                                 remote: $LOOM_REMOTE_BASE/amrex)
#   REMOTE_TMPDIR     TMPDIR for the build                       (default $LOOM_REMOTE_BASE/tmp)
#   EXTRA_MAKE        extra make arguments, word-split           (default none)
#   HIPCC_COMPILE_FLAGS_APPEND, HIPCC_LINK_FLAGS_APPEND  passed into the build
#                     (e.g. the libstdc++ 11 workaround, references/footguns.md)
set -u
usage() { sed -n '2,/^set -u/p' "$0" | sed '$d'; exit "${1:-2}"; }

host=${LOOM_HOST:-local}
# ssh options: LOOM_SSH_OPTS if set (may be empty); otherwise silence the
# weak-crypto warning only where the option exists (OpenSSH >= 10.1; older
# clients reject -o WarnWeakCrypto).
if [ -z "${LOOM_SSH_OPTS+x}" ]; then
    LOOM_SSH_OPTS=""
    read -r ssh_maj ssh_min <<< "$(ssh -V 2>&1 | sed -n 's/^OpenSSH_\([0-9]*\)\.\([0-9]*\).*/\1 \2/p')"
    if [ "${ssh_maj:-0}" -gt 10 ] || { [ "${ssh_maj:-0}" -eq 10 ] && [ "${ssh_min:-0}" -ge 1 ]; }; then
        LOOM_SSH_OPTS="-o WarnWeakCrypto=no"
    fi
fi
# word-split on purpose; bash 3.2 (macOS) rejects "${empty[@]}" under set -u
# shellcheck disable=SC2206
ssh_cmd=(ssh $LOOM_SSH_OPTS "$host")
if [ "$host" = local ]; then
    ssh_cmd=(bash -c)
    B=${LOOM_REMOTE_BASE:-${LOOM_WORK:-$HOME/loom-work}}
    TREE=${REMOTE_TREE:-$(git rev-parse --show-toplevel 2>/dev/null)}
    AMREX_DEFAULT=${AMREX_HOME:-}
else
    B=${LOOM_REMOTE_BASE:-}
    [ -n "$B" ] || { echo "set LOOM_REMOTE_BASE (absolute dir on $host)"; exit 2; }
    TREE=${REMOTE_TREE:-$B/Microphysics}
    AMREX_DEFAULT=$B/amrex
fi
[ -n "$TREE" ] || { echo "set REMOTE_TREE (not in a git tree)"; exit 2; }
SHIP_ROOT=${REMOTE_SHIP_ROOT:-$B/ship}
q() { printf '%q' "$1"; }

apply=0 build=0 ksub=""
case "${1:-}" in
-h|--help) usage 0;;
'') usage;;
status)
    [ $# -eq 1 ] || usage
    "${ssh_cmd[@]}" "L=$(q "$SHIP_ROOT/ship-last-build")
[ -d \"\$L\" ] || { echo \"no build shipped yet (\$L missing)\"; exit 1; }
echo \"ship dir: \$(readlink -f \"\$L\")\"
if [ -e \"\$L/BUILDING\" ]; then
  if kill -0 \"\$(cat \"\$L/BUILDING\")\" 2>/dev/null; then echo '(build still running)'
  else echo '(BUILD DIED without finishing; host up since '\"\$(uptime -s)\"')'; fi
fi
tail -n 25 \"\$L/build.log\" 2>/dev/null"
    exit $?
    ;;
esac
while [ $# -gt 0 ]; do
    case "$1" in
    --apply) apply=1;;
    --build) build=1; apply=1;;
    -h|--help) usage 0;;
    -*) echo "unknown option $1"; usage;;
    *) [ -z "$ksub" ] || usage; ksub=$1;;
    esac
    shift
done
[ -n "$ksub" ] || usage
ksub=${ksub%/}

if [ -n "${LOCAL_TEST_DIR:-}" ]; then
    tdir=$LOCAL_TEST_DIR
else
    top=$(git rev-parse --show-toplevel 2>/dev/null) || { echo "not in a git tree; set LOCAL_TEST_DIR"; exit 2; }
    tdir=$top/unit_test/test_react
fi
tdir=$(cd "$tdir" 2>/dev/null && pwd) || { echo "no local test_react dir: $tdir"; exit 2; }
[ -d "$tdir/$ksub" ] || { echo "KERNEL_SUBDIR $ksub is not a directory in $tdir"; exit 2; }
case "$ksub" in /*|*..*) echo "KERNEL_SUBDIR must be relative to test_react, without .."; exit 2;; esac
if [ $build -eq 1 ]; then
    for v in NETWORK_DIR EOS_DIR INTEGRATOR_DIR; do
        [ -n "${!v:-}" ] || { echo "--build needs $v (no default)"; exit 2; }
    done
    [ -n "${AMREX_HOME:-$AMREX_DEFAULT}" ] || { echo "--build needs AMREX_HOME (AMReX checkout)"; exit 2; }
fi

files=()
for f in GNUmakefile Make.package _parameters main.cpp loom_validate.cpp "$tdir"/loom_*.H; do
    f=${f#"$tdir"/}
    [ -f "$tdir/$f" ] && files+=("$f")
done

stamp=ship-$(date +%Y%m%d-%H%M%S)-$$
stage=$(mktemp -d)
trap 'rm -rf "$stage"' EXIT
mkdir -p "$stage/$stamp/files"
(cd "$tdir" && COPYFILE_DISABLE=1 tar --no-xattrs -cf - "${files[@]}") | (cd "$stage/$stamp/files" && tar -xf -) || exit 1
(cd "$tdir" && COPYFILE_DISABLE=1 tar --no-xattrs --exclude=tests --exclude=build --exclude='build-*' \
    --exclude='tmp_build_dir' --exclude='*.hsaco' --exclude='*.loom' --exclude='.DS_Store' -cf - "$ksub") \
    | (cd "$stage/$stamp/files" && tar -xf -) || exit 1

if [ $build -eq 1 ]; then
    cat > "$stage/$stamp/build.sh" <<EOF
#!/bin/bash
# generated by ship_test_react.sh --build
set -u
S=$(q "$SHIP_ROOT/$stamp")
echo \$\$ > "\$S/BUILDING"
trap 'rm -f "\$S/BUILDING"' EXIT
export ROCM_PATH=$(q "${ROCM_PATH:-/opt/rocm}")
export PATH="\$ROCM_PATH/bin:\$ROCM_PATH/llvm/bin:\$PATH"
export AMREX_HOME=$(q "${AMREX_HOME:-$AMREX_DEFAULT}")
$( [ -n "${HIPCC_COMPILE_FLAGS_APPEND:-}" ] && echo "export HIPCC_COMPILE_FLAGS_APPEND=$(q "$HIPCC_COMPILE_FLAGS_APPEND")" )
$( [ -n "${HIPCC_LINK_FLAGS_APPEND:-}" ] && echo "export HIPCC_LINK_FLAGS_APPEND=$(q "$HIPCC_LINK_FLAGS_APPEND")" )
ARCH=$(q "${ARCH:-}")
[ -n "\$ARCH" ] || ARCH=\$({ amdgpu-arch || rocminfo; } 2>/dev/null | grep -m1 -o 'gfx[0-9a-f]*')
[ -n "\$ARCH" ] || { echo "# WARNING: no ARCH given and no GPU detected; building for gfx90a"; ARCH=gfx90a; }
export MICROPHYSICS_HOME=$(q "$TREE")
export TMPDIR=$(q "${REMOTE_TMPDIR:-$B/tmp}")
mkdir -p "\$TMPDIR"
cd "\$MICROPHYSICS_HOME/unit_test/test_react" || { echo "FAILED: no test_react in \$MICROPHYSICS_HOME"; exit 1; }
echo "# \$(date) host=\$(hostname) dir=\$PWD arch=\$ARCH"
echo "# hipcc: \$(hipcc --version 2>/dev/null | grep -m1 -i 'HIP version')"
echo "# HIPCC_COMPILE_FLAGS_APPEND=\${HIPCC_COMPILE_FLAGS_APPEND:-} HIPCC_LINK_FLAGS_APPEND=\${HIPCC_LINK_FLAGS_APPEND:-}"
t0=\$(date +%s)
nice -n 10 make -j${JOBS:-\$(nproc)} USE_HIP=TRUE AMD_ARCH="\$ARCH" NETWORK_DIR=$(q "$NETWORK_DIR") \\
    EOS_DIR=$(q "$EOS_DIR") INTEGRATOR_DIR=$(q "$INTEGRATOR_DIR") USE_NEUTRINOS=FALSE \\
    XTRA_CXXFLAGS=-ffp-contract=off USE_LOOM_BURN=TRUE LOOM_KERNEL_DIR=$(q "$ksub") ${EXTRA_MAKE:-}
rc=\$?
if [ \$rc -eq 0 ]; then
  echo "# BUILD OK (\$(( \$(date +%s) - t0 ))s): \$(ls -t ./*.ex 2>/dev/null | head -1)"
else
  echo "# BUILD FAILED rc=\$rc (\$(( \$(date +%s) - t0 ))s)"; grep -m10 -E 'error' "\$S/build.log" | cut -c1-300
fi
EOF
fi
(cd "$stage" && COPYFILE_DISABLE=1 tar --no-xattrs -czf "$stage/payload.tgz" "$stamp") || exit 1

"${ssh_cmd[@]}" "set -u
R=$(q "$SHIP_ROOT"); S=\"\$R/$stamp\"; T=$(q "$TREE/unit_test/test_react")
[ -d \"\$T\" ] || { echo \"no test_react in remote tree: \$T\"; exit 1; }
L=\"\$R/ship-last-build\"
if [ $build -eq 1 ] && [ -e \"\$L/BUILDING\" ] && kill -0 \"\$(cat \"\$L/BUILDING\")\" 2>/dev/null; then
  echo \"BUSY: a build is running (\$(readlink -f \"\$L\"))\"; exit 3
fi
mkdir -p \"\$R\" || exit 1
tar xzf - -C \"\$R\" || { echo TAR-FAILED; exit 7; }
cd \"\$S/files\" || exit 8
n_new=0; n_diff=0; n_same=0; : > \"\$S/changed.txt\"
for f in \$(find . -type f | sed 's|^\\./||' | sort); do
  if [ ! -e \"\$T/\$f\" ]; then echo \"NEW     \$f\"; n_new=\$((n_new + 1)); echo \"\$f\" >> \"\$S/changed.txt\"
  elif ! cmp -s \"\$f\" \"\$T/\$f\"; then echo \"DIFFERS \$f\"; n_diff=\$((n_diff + 1)); echo \"\$f\" >> \"\$S/changed.txt\"
  else n_same=\$((n_same + 1)); fi
done
echo \"staged in \$S: \$n_new new, \$n_diff differ, \$n_same identical (vs \$T)\"
if [ $apply -eq 1 ]; then
  while read -r f; do mkdir -p \"\$T/\$(dirname \"\$f\")\" && cp -p \"\$f\" \"\$T/\$f\" || { echo \"COPY FAILED: \$f\"; exit 9; }; done < \"\$S/changed.txt\"
  echo \"applied \$(wc -l < \"\$S/changed.txt\") file(s) to \$T\"
fi
[ $build -eq 1 ] || exit 0
ln -sfn \"\$S\" \"\$L\"
cd \"\$S\" || exit 8
nohup setsid bash build.sh > build.log 2>&1 < /dev/null &
sleep 1; echo \"build launched at \$(date); log: \$S/build.log (ship_test_react.sh status)\"" < "$stage/payload.tgz"

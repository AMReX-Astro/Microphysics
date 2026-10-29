#!/usr/bin/env python3
"""Generate a Loom kernel's per-network data headers from the stock code.

For a network that uses the networks/rhs.H framework (NEW_NETWORK_IMPLEMENTATION;
aprox13, aprox19, aprox21, ignition_chamulak, ignition_simple, iso7, powerlaw,
rprox, triple_alpha_plus_cago) this compiles extract_network_data.cpp against
the stock headers with one build's flags, runs it, and writes

  network_data.H  namespace mpl: NumSpec/NumAux, Species::/Rates:: indices,
                  the zion/aion/aion_inv arrays, bion/mion, the rhs_data fields as
                  constexpr if-chain functions (rd_*), the intermediate-rate
                  index table, the screen factors of every pair the network's
                  apply_screening folds, constants the glue uses, and
                  MPL_GEN_* configuration macros.
  burn.H          neqs, net_ienuc, S* indices, JacNetArray2D and NetYArray1D
                  (rhs_state_t::y) on support.H's Array2D/Array1D, and the
                  burn_t mirror (same fields, order and defaults as
                  interfaces/burn_type.H for this configuration).

Both are Loom-compatible (README rules 6/7: constexpr scalar functions with
if-chains; namespace-scope arrays only as internal-linkage constexpr tables)
and every floating-point literal is the
shortest decimal that round-trips to the stock bits.

usage (from anywhere; paths are resolved against the Microphysics root):
  uv run --no-project gen_network_data.py --out DIR [--emit-dir DIR] \
      [--flags FLAGS.txt] [--test-dir unit_test/test_react] \
      [--support-dir DIR] [--cxx clang++] [-- MAKEVAR=VALUE ...]

Flags come from --flags (an inventory.py flags.txt; its relative -I paths are
relative to --test-dir) or, without --flags, from `make -n` with the given
make variables, after generating the build's headers into DIR/build
(TMP_BUILD_DIR, so no other build directory is touched).

Checks (any failure exits non-zero, nothing is left half-written in
--emit-dir): the extractor's own static_asserts; every enumerator value; a
host TU compiled with the stock flags that static_asserts each generated value
against the stock constexpr source (rhs_data fields, intermediate indices,
screen factors, zion/aion/aion_inv, bion/mion, sizes, S* indices, burn_t size
and every field offset); and a standalone `-std=c++20 -fsyntax-only` of the
generated headers with support.H only.
"""
import argparse
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
def _microphysics_root():
    """$MICROPHYSICS_HOME, else the repo this skill sits in
    (<repo>/.claude/skills/<skill> or <repo>/.agents/skills/<skill>), else the
    git checkout of the current directory (a user-level skill install)."""
    def ok(d):
        return d and os.path.isdir(os.path.join(d, 'unit_test', 'test_react'))
    env = os.environ.get('MICROPHYSICS_HOME')
    if ok(env):
        return os.path.abspath(env)
    here = os.path.abspath(os.path.join(SKILL, '..', '..', '..'))
    if ok(here):
        return here
    try:
        top = subprocess.run(['git', 'rev-parse', '--show-toplevel'], capture_output=True,
                             text=True).stdout.strip()
    except OSError:
        top = ''
    if ok(top):
        return top
    sys.exit('cannot find the Microphysics tree: set MICROPHYSICS_HOME or run from inside it')


ROOT = _microphysics_root()


def die(msg):
    sys.exit(f'gen_network_data: {msg}')


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


# ---------------------------------------------------------------------------
# flags
# ---------------------------------------------------------------------------

def flags_from_make(test, build, makevars, out):
    mk = ['make', f'TMP_BUILD_DIR={build}'] + [v for v in makevars
                                               if not v.startswith('XTRA_CXXFLAGS=')]
    dry = run(mk + ['-n'], cwd=test)
    line = next((l for l in dry.stdout.splitlines()
                 if re.search(r' -c main\.cpp( |$)', l)), None)
    if line is None:
        open(os.path.join(out, 'make-n.log'), 'w').write(dry.stdout + dry.stderr)
        die(f'no main.cpp compile line in `make -n` output (see {out}/make-n.log)')
    toks = shlex.split(line)
    obj = toks[toks.index('-o') + 1]
    # generate network_properties.H / extern_parameters.H / AMReX_Config.H by
    # building main.o syntax-only (as inventory.py does); a fresh
    # TMP_BUILD_DIR can race on the generated sources, so retry once
    for _ in range(2):
        gen = run(mk + ['XTRA_CXXFLAGS=-fsyntax-only', obj], cwd=test)
        open(os.path.join(out, 'make.log'), 'a').write(gen.stdout + gen.stderr)
        if any(f == 'network_properties.H' for _, _, fs in os.walk(build) for f in fs):
            break
    else:
        die(f'header generation failed, see {out}/make.log')
    flags, i = [], 1
    while i < len(toks):
        t = toks[i]
        if t.startswith(('-D', '-U', '-I', '-std=')):
            flags.append(t)
        elif t in ('-include', '-isystem', '-iquote'):
            flags += toks[i:i + 2]
            i += 1
        i += 1
    open(os.path.join(out, 'flags.txt'), 'w').write('\n'.join(flags) + '\n')
    return flags


def include_dirs(flags, test):
    dirs, i = [], 0
    while i < len(flags):
        f = flags[i]
        if f.startswith('-I') and len(f) > 2:
            dirs.append(f[2:])
        elif f in ('-isystem', '-iquote') and i + 1 < len(flags):
            dirs.append(flags[i + 1])
            i += 1
        i += 1
    return [d if os.path.isabs(d) else os.path.normpath(os.path.join(test, d)) for d in dirs]


def find_header(dirs, name, prefer=None):
    hits = [os.path.join(d, name) for d in dirs if os.path.isfile(os.path.join(d, name))]
    if prefer:
        pref = [h for h in hits if prefer(h)]
        if pref:
            return pref[0]
    return hits[0] if hits else None


# ---------------------------------------------------------------------------
# enumerator names (values come from the extractor)
# ---------------------------------------------------------------------------

def strip_comments(t):
    t = re.sub(r'/\*.*?\*/', '', t, flags=re.S)
    return re.sub(r'//[^\n]*', '', t)


def enum_names(text, enum):
    m = re.search(r'enum\s+' + enum + r'\b[^{]*\{(.*?)\}', strip_comments(text), re.S)
    if not m:
        return None
    names = []
    for item in m.group(1).split(','):
        item = item.strip()
        if not item:
            continue
        names.append(re.match(r'([A-Za-z_]\w*)', item).group(1))
    return names


# ---------------------------------------------------------------------------
# extractor output
# ---------------------------------------------------------------------------

def hexf(s):
    return float.fromhex(s)


def parse_extract(text):
    d = {'def': {}, 'cfg': {}, 'species': [], 'rates': [], 'ion': {}, 'bm': {},
         'const': {}, 'rd': {}, 'inter': {}, 'densexp': {}, 'screen': [],
         'screen_method': None}
    end = False
    for line in text.splitlines():
        p = line.split()
        if not p:
            continue
        tag = p[0]
        if tag == 'NETWORK':
            d['network'] = p[1]
        elif tag == 'DEF':
            d['def'][p[1]] = int(p[2])
        elif tag == 'SCREEN_METHOD':
            d['screen_method'] = p[1]
        elif tag == 'CFG':
            d['cfg'][p[1]] = int(p[2])
        elif tag == 'SPECIES':
            d['species'].append((p[1], int(p[2])))
        elif tag == 'RATE':
            d['rates'].append((p[1], int(p[2])))
        elif tag == 'ION':
            d['ion'][int(p[1])] = tuple(hexf(x) for x in p[2:7])
        elif tag == 'BIONMION':
            d['bm'][int(p[1])] = (hexf(p[2]), hexf(p[3]))
        elif tag == 'CONST':
            d['const'][p[1]] = hexf(p[2])
        elif tag == 'RD':
            rate, field, vals = int(p[1]), p[2], p[3:]
            conv = hexf if field in ('number', 'forward_branching_ratio',
                                     'reverse_branching_ratio') else int
            v = [conv(x) for x in vals]
            d['rd'].setdefault(rate, {})[field] = v if len(v) > 1 else v[0]
        elif tag == 'INTER':
            d['inter'][int(p[1])] = int(p[2])
        elif tag == 'DENSEXP':
            d['densexp'][int(p[1])] = (int(p[2]), int(p[3]))
        elif tag == 'SCREEN':
            rate, k = int(p[1]), int(p[2])
            z1, a1, z2, a2 = (hexf(x) for x in p[3:7])
            fields = [('z1', z1), ('z2', z2), ('a1', a1), ('a2', a2)]
            for kv in p[7:]:
                n, v = kv.split('=')
                fields.append((n, hexf(v)))
            d['screen'].append((rate, k, (z1, a1, z2, a2), fields))
        elif tag == 'END':
            end = True
    if not end:
        die('extractor output truncated (no END record)')
    return d


# ---------------------------------------------------------------------------
# formatting
# ---------------------------------------------------------------------------

def lit(x):
    """Shortest decimal double literal that round-trips to the same bits."""
    if math.isnan(x) or math.isinf(x):
        die(f'non-finite value {x!r} cannot be written as a literal')
    s = repr(float(x))
    if not any(c in s for c in '.eE'):
        s += '.0'
    back = float(s)
    if back != x or math.copysign(1.0, back) != math.copysign(1.0, x):
        die(f'literal {s} does not round-trip')
    return s


def int_of(x, what):
    if x != int(x):
        die(f'{what} = {x!r} is not integral')
    return int(x)


RD_DEFAULTS = {  # interfaces/rhs_type.H rhs_t member defaults
    'species': [-1] * 6, 'number': [0.0] * 6, 'exponent': [0] * 6,
    'forward_branching_ratio': 1.0, 'reverse_branching_ratio': 1.0,
    'apply_identical_particle_factor': 1, 'rate_can_be_tabulated': 1,
    'screen_forward_reaction': 1, 'screen_reverse_reaction': 1,
    'additional_reaction': [-1] * 3,
}

SCREEN_FIELDS = {  # screening/screen_data.H screen_factors_t, per SCREEN_METHOD
    'screen5': ['zs53', 'z1_53', 'z2_53', 'zs13', 'zs13inv', 'zhat', 'zhat2', 'lzav', 'aznut'],
    'chugunov2007': ['ztilde', 'ztilde3'],
    'chugunov2009': ['zs52', 'z1_52', 'z2_52', 'zs53', 'z1_53', 'z2_53', 'aznut', 'ztilde'],
    'chabrier1998': ['zs53', 'z1_53', 'z2_53', 'zs13', 'zs13inv', 'aznut'],
}


class Gen:
    def __init__(self, d, cmdline, flags_src):
        self.d, self.cmdline, self.flags_src = d, cmdline, flags_src
        cfg = d['cfg']
        self.nspec, self.ntot, self.naux = cfg['NumSpec'], cfg['NumSpecTotal'], cfg['NumAux']
        self.nrates = cfg['NumRates']
        self.spname = {v: n for n, v in d['species']}
        self.rname = {v: n for n, v in d['rates']}

    # names for indices in emitted code
    def sp(self, v):
        return 'Species::' + self.spname[v] if v in self.spname else str(v)

    def rt(self, v):
        return 'Rates::' + self.rname[v] if v in self.rname else str(v)

    def provenance(self, what):
        d = self.d
        defs = ' '.join(f'{k}={v}' for k, v in sorted(d['def'].items()))
        return [
            f'// {what} for the Loom burn kernel: network {d["network"]}.',
            '// GENERATED by the loom-burn-kernel skill, scripts/gen_network_data.py',
            '// from the stock headers (extract_network_data.cpp). Do not edit; rerun:',
            f'//   {self.cmdline}',
            f'// flags: {self.flags_src}',
            f'// defines: {defs}' + (f' SCREEN_METHOD={d["screen_method"]}' if d['screen_method'] else ''),
        ]

    def screen_pairs(self):
        """Unique (Z1, A1, Z2, A2) pairs in first-use order, with their fields."""
        pairs, seen = [], {}
        for rate, k, key, fields in self.d['screen']:
            if key in seen:
                if pairs[seen[key]][1] != fields:
                    die(f'screen pair {key} has different factors for rate {rate}')
                continue
            seen[key] = len(pairs)
            pairs.append((key, fields, (rate, k)))
        return pairs

    def network_data(self):
        d, L = self.d, []
        guard = 'MPL_NETWORK_DATA_H'
        L += self.provenance('Per-network data (networks/rhs.H framework)')
        L += ['//',
              '// Generic code (network_rhs.H, screen.H, eos.H, burn.H) reads everything',
              '// network-specific through this header; rhs_data(rate) is split into one',
              '// constexpr if-chain function per rhs_t field (README rule 6); aion[],',
              '// zion[], aion_inv[] are namespace-scope constexpr arrays as in the stock',
              '// HIP build (rule 7).',
              f'#ifndef {guard}', f'#define {guard}', '', '#include "support.H"', '']
        sm = d['screen_method'] or 'none'
        L += [f'#define MPL_GEN_NETWORK_NAME "{d["network"]}"',
              f'#define MPL_GEN_SCREENING {d["def"]["SCREENING"]}',
              f'#define MPL_GEN_SCREEN_METHOD_{sm} 1',
              f'#define MPL_GEN_NEUTRINOS {d["def"]["NEUTRINOS"]}',
              f'#define MPL_GEN_NSE {d["def"]["NSE"]}',
              f'#define MPL_GEN_NSE_NET {d["def"]["NSE_NET"]}',
              f'#define MPL_GEN_NSE_TABLE {d["def"]["NSE_TABLE"]}',
              f'#define MPL_GEN_STRANG {d["def"]["STRANG"]}',
              f'#define MPL_GEN_SIMPLIFIED_SDC {d["def"]["SIMPLIFIED_SDC"]}',
              f'#define MPL_GEN_TRUE_SDC {d["def"]["TRUE_SDC"]}',
              f'#define MPL_GEN_HAVE_BION_MION {d["def"]["HAVE_BION_MION"]}',
              '']
        ion = d['ion']
        def table(name, k):
            return [f'constexpr Real {name}[NumSpec + NumSpecExtra] = {{',
                    '    ' + ', '.join(lit(ion[s][k]) for s in range(1, self.ntot + 1)), '};']
        L += ['namespace mpl {', '',
              '// network_properties.H',
              f'constexpr int NumSpec = {self.nspec};',
              f'constexpr int NumSpecExtra = {d["cfg"]["NumSpecExtra"]};',
              f'constexpr int NumAux = {self.naux};', '',
              '// network_properties.H aion[] / aion_inv[] / zion[] (0-based, NumSpecTotal',
              '// entries): HIP_CONSTEXPR arrays in an anonymous namespace, as in the stock',
              '// HIP build (namespace-scope constexpr arrays import, probes d01, d16)',
              'namespace {']
        L += table('aion', 2) + table('aion_inv', 4) + table('zion', 3)
        L += ['} // namespace', '',
              '// network_properties.H namespace Species (NumSpecTotal = last species)',
              'namespace Species {']
        L += [f'constexpr int {n} = {v};' for n, v in d['species']]
        L += [f'constexpr int NumSpecTotal = {self.ntot};', '} // namespace Species', '',
              '// actual_network.H Rates::NetworkRates',
              'namespace Rates {']
        L += [f'constexpr int {n} = {v};' for n, v in d['rates']]
        L += [f'constexpr int NumRates = {self.nrates};', '} // namespace Rates', '']

        def chain(sig, items, default, comment=None):
            out = ([comment] if comment else []) + [f'MPL_DEVICE constexpr {sig}', '{']
            out += [f'    if ({c}) {{ return {v}; }}' for c, v in items]
            out += [f'    return {default};', '}', '']
            return out

        sps = range(1, self.ntot + 1)
        L += chain('Real net_aion(int spec)', [(f'spec == {self.sp(s)}', lit(ion[s][0])) for s in sps], '-1.0',
                   '// NetworkProperties::aion(spec) / zion(spec) (-1 for an unknown species)')
        L += chain('Real net_zion(int spec)', [(f'spec == {self.sp(s)}', lit(ion[s][1])) for s in sps], '-1.0')
        L += chain('Real net_aion_inv(int spec)', [(f'spec == {self.sp(s)}', lit(ion[s][4])) for s in sps], '-1.0',
                   '// network_properties.H aion_inv[spec - 1]')
        integral = all(v == int(v) for s in sps for v in ion[s][:2])
        if integral:
            L += chain('int net_aion_i(int spec)', [(f'spec == {self.sp(s)}', str(int(ion[s][0]))) for s in sps],
                       '-1', '// integer A and Z (all integral in this network), usable as template arguments')
            L += chain('int net_zion_i(int spec)', [(f'spec == {self.sp(s)}', str(int(ion[s][1]))) for s in sps],
                       '-1')
        else:
            L += ['// net_aion_i / net_zion_i: not generated, A or Z is not integral for some species', '']
        if d['def']['HAVE_BION_MION']:
            bm = d['bm']
            ns = range(1, self.nspec + 1)
            L += chain('Real net_bion(int spec)', [(f'spec == {self.sp(s)}', lit(bm[s][0])) for s in ns], '0.0',
                       '// actual_network.H network::bion<spec>() / mion<spec>(), spec = 1..NumSpec')
            L += chain('Real net_mion(int spec)', [(f'spec == {self.sp(s)}', lit(bm[s][1])) for s in ns], '0.0')
        else:
            L += ['// net_bion / net_mion: this network defines no network::bion<>() / mion<>()', '']
        c = d['const']
        L += ['// fundamental_constants.H values the per-network glue (ener_gener_rate) uses',
              f'constexpr Real net_legacy_enuc_conv2 = {lit(c["legacy_enuc_conv2"])};   // C::Legacy::enuc_conv2',
              f'constexpr Real net_legacy_n_A = {lit(c["legacy_n_A"])};   // C::Legacy::n_A',
              f'constexpr Real net_legacy_MeV2erg = {lit(c["legacy_MeV2erg"])};   // C::Legacy::MeV2erg',
              f'constexpr Real net_legacy_MeV2gr = {lit(c["legacy_MeV2gr"])};   // C::Legacy::MeV2gr',
              f'constexpr Real net_enuc_conv2 = {lit(c["enuc_conv2"])};   // C::enuc_conv2',
              f'constexpr Real net_n_A = {lit(c["n_A"])};   // C::n_A', '']

        # rhs_data
        L += ['// ---------------------------------------------------------------------------',
              '// actual_network.H RHS::rhs_data(rate), one function per rhs_t field',
              '// (interfaces/rhs_type.H defaults for rates not listed). slot 0..5 = A..F.',
              '// ---------------------------------------------------------------------------', '',
              'MPL_DEVICE constexpr int rd_pick_i(int slot, int a, int b, int c, int d, int e, int f)',
              '{', '    if (slot == 0) { return a; }', '    if (slot == 1) { return b; }',
              '    if (slot == 2) { return c; }', '    if (slot == 3) { return d; }',
              '    if (slot == 4) { return e; }', '    return f;', '}', '',
              'MPL_DEVICE constexpr Real rd_pick_r(int slot, Real a, Real b, Real c, Real d, Real e, Real f)',
              '{', '    if (slot == 0) { return a; }', '    if (slot == 1) { return b; }',
              '    if (slot == 2) { return c; }', '    if (slot == 3) { return d; }',
              '    if (slot == 4) { return e; }', '    return f;', '}', '',
              'MPL_DEVICE constexpr int rd_pick3(int k, int a, int b, int c)',
              '{', '    if (k == 1) { return a; }', '    if (k == 2) { return b; }', '    return c;', '}', '']
        rd = d['rd']
        rates = range(1, self.nrates + 1)

        def rows(field, fmt):
            return [(f'rate == {self.rt(r)}', fmt(rd[r][field])) for r in rates
                    if rd[r][field] != RD_DEFAULTS[field]]

        L += chain('int rd_species(int rate, int slot)',
                   rows('species', lambda v: 'rd_pick_i(slot, ' + ', '.join(self.sp(x) if x >= 0 else str(x) for x in v) + ')'),
                   '-1')
        L += chain('Real rd_number(int rate, int slot)',
                   rows('number', lambda v: 'rd_pick_r(slot, ' + ', '.join(lit(x) for x in v) + ')'), '0.0')
        L += chain('int rd_exponent(int rate, int slot)',
                   rows('exponent', lambda v: 'rd_pick_i(slot, ' + ', '.join(str(x) for x in v) + ')'), '0')
        L += chain('Real rd_forward_branching_ratio(int rate)', rows('forward_branching_ratio', lit), '1.0')
        L += chain('Real rd_reverse_branching_ratio(int rate)', rows('reverse_branching_ratio', lit), '1.0')
        for f in ('apply_identical_particle_factor', 'rate_can_be_tabulated',
                  'screen_forward_reaction', 'screen_reverse_reaction'):
            L += chain(f'int rd_{f}(int rate)', rows(f, str), str(RD_DEFAULTS[f]))
        L += chain('int rd_additional_reaction(int rate, int k)',
                   rows('additional_reaction',
                        lambda v: 'rd_pick3(k, ' + ', '.join(self.rt(x) if x >= 1 else str(x) for x in v) + ')'),
                   '-1')

        # intermediate rates
        inter = d['inter']
        L += ['// networks/rhs.H num_intermediate_reactions() and',
              '// locate_intermediate_rate_index(rate) (-1: not intermediate)',
              f'constexpr int net_num_intermediate = {d["cfg"]["num_intermediate"]};', '']
        L += chain('int net_intermediate_index(int rate)',
                   [(f'rate == {self.rt(r)}', str(inter[r])) for r in rates if inter[r] >= 1], '-1')

        # screening
        if d['def']['SCREENING']:
            method = d['screen_method']
            if method not in SCREEN_FIELDS:
                die(f'SCREEN_METHOD {method} has no field list in gen_network_data.py')
            pairs = self.screen_pairs()
            fields = ['z1', 'z2', 'a1', 'a2'] + SCREEN_FIELDS[method]
            L += ['// ---------------------------------------------------------------------------',
                  '// screening/screen_data.H calculate_screen_factor(Z1, A1, Z2, A2) as networks/',
                  '// rhs.H apply_screening folds it (constexpr) for every screened pair of this',
                  f'// network ({method} fields). net_screen_pair(Z1, A1, Z2, A2) -> pair index',
                  '// (-1: not screened by this network); net_screen_<field>(pair) -> value.',
                  '// ---------------------------------------------------------------------------',
                  f'constexpr int net_num_screen_pairs = {len(pairs)};', '']
            out = ['MPL_DEVICE constexpr int net_screen_pair(int Z1, int A1, int Z2, int A2)', '{']
            for i, (key, _, (rate, k)) in enumerate(pairs):
                z1, a1, z2, a2 = (int_of(v, f'screened Z/A of rate {self.rname.get(rate, rate)}') for v in key)
                note = self.rname.get(rate, str(rate)) + (', second factor' if k == 2 else '')
                out.append(f'    if (Z1 == {z1} && A1 == {a1} && Z2 == {z2} && A2 == {a2}) {{ return {i}; }}  // {note}')
            out += ['    return -1;', '}', '']
            L += out
            for f in fields:
                L += chain(f'Real net_screen_{f}(int pair)',
                           [(f'pair == {i}', lit(dict(fl)[f])) for i, (_, fl, _) in enumerate(pairs)], '0.0')
        L += ['} // namespace mpl', '', '#endif', '']
        return '\n'.join(L)

    def burn_fields(self):
        """(type, name, array_len or None, default) in interfaces/burn_type.H order."""
        d = self.d['def']
        R = 'Real'
        f = [(R, 'time', None, ''), (R, 'rho', None, ''), (R, 'T', None, ''), (R, 'e', None, ''),
             (R, 'xn', 'NumSpec', '')]
        if self.naux > 0:
            f.append((R, 'aux', 'NumAux', ''))
        f += [(R, 'e_scale', None, ''), (R, 'dedr', None, ''), (R, 'dedT', None, ''),
              (R, 'dedA', None, ''), (R, 'dedZ', None, ''),
              ('int', 'i', None, ''), ('int', 'j', None, ''), ('int', 'k', None, ''),
              (R, 'y', 'SVAR', ''), (R, 'rho_orig', None, ''), (R, 'rhoe_orig', None, ''),
              (R, 'ydot_a', 'SVAR', ''), ('int', 'sdc_iter', None, ''), ('int', 'num_sdc_iters', None, ''),
              (R, 'T_fixed', None, '-1.0'), (R, 'cv', None, ''), (R, 'dx', None, ''),
              (R, 'mu', None, ''), (R, 'mu_e', None, ''), (R, 'y_e', None, ''), (R, 'eta', None, ''),
              (R, 'abar', None, ''), (R, 'zbar', None, '')]
        if d['NSE_NET']:
            f += [(R, 'mu_p', None, ''), (R, 'mu_n', None, '')]
        if d['NSE']:
            f.append(('bool', 'nse', None, ''))
        f += [('int', 'n_rhs', None, ''), ('int', 'n_jac', None, ''), ('int', 'n_step', None, ''),
              ('bool', 'success', None, ''), ('short', 'error_code', None, '')]
        return f

    def burn(self):
        d, cfg = self.d['def'], self.d['cfg']
        L = self.provenance('Shared burn state')
        conf = [f'NAUX_NET={self.naux}', 'STRANG' if d['STRANG'] else
                ('SIMPLIFIED_SDC' if d['SIMPLIFIED_SDC'] else ('TRUE_SDC' if d['TRUE_SDC'] else 'no coupling define')),
                'NSE' if d['NSE'] else 'no NSE', 'NSE_NET' if d['NSE_NET'] else 'no NSE_NET', 'no NONAKA_PLOT']
        L += ['//',
              '// from interfaces/burn_type.H burn_t, for this configuration: ' + ', '.join(conf) + '.',
              '// Same fields, same order, same default values as the stock struct;',
              '// gen_network_data.py static_asserts size and every offset against it.',
              '#ifndef MPL_BURN_H', '#define MPL_BURN_H', '', '#include "support.H"',
              '#include "network_data.H"', '', 'namespace mpl {', '',
              '// interfaces/burn_type.H']
        for n in ('neqs', 'net_ienuc', 'SFS', 'SEINT', 'SFX', 'SRHO', 'SMX', 'SMY', 'SMZ',
                  'SEDEN', 'SVAR', 'SVAR_EVOLVE'):
            L.append(f'constexpr int {n} = {cfg[n]};')
        L += ['', '// interfaces/burn_type.H JacNetArray2D (jac_t = Real) and rhs_type.H rhs_state_t::y',
              'using JacNetArray2D = Array2D<Real, 1, neqs, 1, neqs>;',
              'using NetYArray1D = Array1D<Real, 1, NumSpec>;']
        L += ['', 'struct burn_t', '{']
        for t, n, arr, dflt in self.burn_fields():
            L.append(f'    {t} {n}' + (f'[{arr}]' if arr else '') + '{' + dflt + '};')
        L += ['};', '', '} // namespace mpl', '', '#endif', '']
        return '\n'.join(L)

    def check_tu(self, emit):
        """Host TU (stock flags): static_assert every generated value against the stock source."""
        d = self.d
        L = ['// generated by gen_network_data.py: generated headers vs stock constexpr source',
             '#include <cstddef>', '#include <AMReX_Loop.H>', '#include <network.H>',
             '#include <burn_type.H>', '#include <fundamental_constants.H>', '#include <rhs.H>',
             '#ifdef SCREENING', '#include <screen_data.H>', '#endif',
             f'#include "{os.path.join(emit, "network_data.H")}"',
             f'#include "{os.path.join(emit, "burn.H")}"', '',
             'static_assert(NumSpec == mpl::NumSpec && NumAux == mpl::NumAux && NumSpecExtra == mpl::NumSpecExtra);',
             'static_assert(NumSpecTotal == mpl::Species::NumSpecTotal && Rates::NumRates == mpl::Rates::NumRates);',
             'static_assert(RHS::num_intermediate_reactions() == mpl::net_num_intermediate);']
        for n in ('neqs', 'net_ienuc', 'SFS', 'SEINT', 'SFX', 'SRHO', 'SMX', 'SMY', 'SMZ',
                  'SEDEN', 'SVAR', 'SVAR_EVOLVE'):
            L.append(f'static_assert({n} == mpl::{n});')
        L += [f'static_assert(Species::{n} == mpl::Species::{n});' for n, _ in d['species']]
        L += [f'static_assert(Rates::{n} == mpl::Rates::{n});' for n, _ in d['rates']]
        L += ['static_assert(sizeof(burn_t) == sizeof(mpl::burn_t), "burn_t size");']
        L += [f'static_assert(offsetof(burn_t, {n}) == offsetof(mpl::burn_t, {n}), "burn_t.{n}");'
              for _, n, _, _ in self.burn_fields()]
        L += ['static_assert(C::Legacy::enuc_conv2 == mpl::net_legacy_enuc_conv2 && C::Legacy::n_A == mpl::net_legacy_n_A);',
              'static_assert(C::Legacy::MeV2erg == mpl::net_legacy_MeV2erg && C::Legacy::MeV2gr == mpl::net_legacy_MeV2gr);',
              'static_assert(C::enuc_conv2 == mpl::net_enuc_conv2 && C::n_A == mpl::net_n_A);', '',
              'template <int rate>', 'constexpr bool rd_ok()', '{',
              '    constexpr RHS::rhs_t d = RHS::rhs_data(rate);',
              '    return d.species_A == mpl::rd_species(rate, 0) && d.species_B == mpl::rd_species(rate, 1) &&',
              '           d.species_C == mpl::rd_species(rate, 2) && d.species_D == mpl::rd_species(rate, 3) &&',
              '           d.species_E == mpl::rd_species(rate, 4) && d.species_F == mpl::rd_species(rate, 5) &&',
              '           d.number_A == mpl::rd_number(rate, 0) && d.number_B == mpl::rd_number(rate, 1) &&',
              '           d.number_C == mpl::rd_number(rate, 2) && d.number_D == mpl::rd_number(rate, 3) &&',
              '           d.number_E == mpl::rd_number(rate, 4) && d.number_F == mpl::rd_number(rate, 5) &&',
              '           d.exponent_A == mpl::rd_exponent(rate, 0) && d.exponent_B == mpl::rd_exponent(rate, 1) &&',
              '           d.exponent_C == mpl::rd_exponent(rate, 2) && d.exponent_D == mpl::rd_exponent(rate, 3) &&',
              '           d.exponent_E == mpl::rd_exponent(rate, 4) && d.exponent_F == mpl::rd_exponent(rate, 5) &&',
              '           d.forward_branching_ratio == mpl::rd_forward_branching_ratio(rate) &&',
              '           d.reverse_branching_ratio == mpl::rd_reverse_branching_ratio(rate) &&',
              '           d.apply_identical_particle_factor == mpl::rd_apply_identical_particle_factor(rate) &&',
              '           d.rate_can_be_tabulated == mpl::rd_rate_can_be_tabulated(rate) &&',
              '           d.screen_forward_reaction == mpl::rd_screen_forward_reaction(rate) &&',
              '           d.screen_reverse_reaction == mpl::rd_screen_reverse_reaction(rate) &&',
              '           d.additional_reaction_1 == mpl::rd_additional_reaction(rate, 1) &&',
              '           d.additional_reaction_2 == mpl::rd_additional_reaction(rate, 2) &&',
              '           d.additional_reaction_3 == mpl::rd_additional_reaction(rate, 3) &&',
              '           RHS::locate_intermediate_rate_index(rate) == mpl::net_intermediate_index(rate);',
              '}', '']
        L += [f'static_assert(rd_ok<{r}>(), "rhs_data / intermediate index of rate {r}");'
              for r in range(1, self.nrates + 1)]
        for s in range(1, self.ntot + 1):
            L.append(f'static_assert(NetworkProperties::aion({s}) == mpl::net_aion({s}) && '
                     f'NetworkProperties::zion({s}) == mpl::net_zion({s}));')
        if d['def']['HAVE_BION_MION']:
            for s in range(1, self.nspec + 1):
                L.append(f'static_assert(network::bion<{s}>() == mpl::net_bion({s}) && '
                         f'network::mion<{s}>() == mpl::net_mion({s}));')
        L += ['// aion[] / zion[] / aion_inv[] are runtime arrays on the CPU build: checked in main()',
              'int main()', '{', '    int bad = 0;',
              '    for (int n = 0; n < NumSpecTotal; ++n) {',
              '        bad += (mpl::aion[n] != aion[n]) + (mpl::zion[n] != zion[n]) +',
              '               (mpl::aion_inv[n] != aion_inv[n]);',
              '    }',
              '    for (int s = 1; s <= NumSpecTotal; ++s) {',
              '        bad += (mpl::net_aion_inv(s) != aion_inv[s - 1]);',
              '    }']
        if d['def']['SCREENING']:
            fields = ['z1', 'z2', 'a1', 'a2'] + SCREEN_FIELDS[d['screen_method']]
            for i, (key, _, _) in enumerate(self.screen_pairs()):
                args = ', '.join(lit(v) for v in key)
                cond = ' && '.join(f'sf{i}.{f} == mpl::net_screen_{f}({i})' for f in fields)
                L.insert(L.index('// aion[] / zion[] / aion_inv[] are runtime arrays on the CPU build: checked in main()'),
                         f'constexpr auto sf{i} = scrn::calculate_screen_factor({args});\n'
                         f'static_assert({cond}, "screen pair {i}");\n'
                         f'static_assert(mpl::net_screen_pair({", ".join(str(int(v)) for v in key)}) == {i});')
        L += ['    return bad;', '}', '']
        return '\n'.join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--out', required=True, help='work directory (extractor, logs, flags.txt)')
    ap.add_argument('--emit-dir', help='where network_data.H and burn.H go (default: --out)')
    ap.add_argument('--flags', help='flags.txt (from inventory.py or a previous run)')
    ap.add_argument('--test-dir', default='unit_test/test_react')
    ap.add_argument('--support-dir', default=os.path.join(SKILL, 'assets', 'kernel'),
                    help='directory with support.H for the standalone syntax check')
    ap.add_argument('--cxx', default='clang++')
    ap.add_argument('makevars', nargs='*')
    a = ap.parse_args()

    out = os.path.abspath(a.out)
    emit = os.path.abspath(a.emit_dir or a.out)
    a.support_dir = os.path.abspath(a.support_dir)
    test = a.test_dir if os.path.isabs(a.test_dir) else os.path.join(ROOT, a.test_dir)
    os.makedirs(out, exist_ok=True)
    os.makedirs(emit, exist_ok=True)
    cmdline = 'uv run --no-project ' + ' '.join(shlex.quote(x) for x in
                                                [os.path.relpath(__file__, ROOT)] + sys.argv[1:])

    if a.flags:
        flags = open(a.flags).read().split()
        flags_src = os.path.relpath(os.path.abspath(a.flags), ROOT)
    else:
        if not a.makevars:
            die('give --flags FILE or make variables after `--` (e.g. NETWORK_DIR=iso7 ...)')
        flags = flags_from_make(test, os.path.join(out, 'build'), a.makevars, out)
        flags_src = 'make -n ' + ' '.join(a.makevars)

    dirs = include_dirs(flags, test)
    net_h = find_header(dirs, 'actual_network.H', lambda h: '/networks/' in h)
    props = find_header(dirs, 'network_properties.H')
    if not net_h or not props:
        die('actual_network.H / network_properties.H not on the include path '
            '(were the build headers generated?)')
    net_text = open(net_h).read()
    if 'NEW_NETWORK_IMPLEMENTATION' not in net_text or not re.search(r'\brhs_data\s*\(', net_text):
        die(f'{net_h} does not use the networks/rhs.H framework '
            '(no NEW_NETWORK_IMPLEMENTATION / RHS::rhs_data); not supported')
    rates = enum_names(net_text, 'NetworkRates')
    species = enum_names(open(props).read(), 'NetworkSpecies')
    if rates is None or species is None:
        die('could not find enum NetworkRates / NetworkSpecies')
    rates = [n for n in rates if n != 'NumRates']
    species = [n for n in species if n not in ('NumberSpecies', 'NumberExtraSpecies', 'NumberTotalSpecies')]
    have_bm = int(bool(re.search(r'constexpr\s+amrex::Real\s+bion\s*\(\s*\)', net_text)) and
                  bool(re.search(r'constexpr\s+amrex::Real\s+mion\s*\(\s*\)', net_text)))

    with open(os.path.join(out, 'extract_names.inc'), 'w') as f:
        f.writelines(f'MPL_X_SPECIES({n})\n' for n in species)
        f.writelines(f'MPL_X_RATE({n})\n' for n in rates)
    exe = os.path.join(out, 'extract_network_data')
    cc = [a.cxx, *flags, f'-I{out}', f'-DMPL_X_HAVE_BION_MION={have_bm}', '-O0', '-w',
          os.path.join(HERE, 'extract_network_data.cpp'), '-o', exe]
    r = run(cc, cwd=test)
    open(os.path.join(out, 'extract_build.log'), 'w').write(' '.join(map(shlex.quote, cc)) + '\n' +
                                                           r.stdout + r.stderr)
    if r.returncode:
        die(f'extractor failed to build, see {out}/extract_build.log:\n' +
            '\n'.join(l for l in r.stderr.splitlines() if 'error' in l)[:3000])
    r = run([exe], cwd=test)
    open(os.path.join(out, 'extract.txt'), 'w').write(r.stdout)
    if r.returncode:
        die(f'extractor exited {r.returncode}: {r.stderr[-2000:]}')
    d = parse_extract(r.stdout)

    # consistency of what was extracted, and what the kernel supports
    cfg, df = d['cfg'], d['def']
    if cfg['sizeof_Real'] != 8:
        die('amrex::Real is not double (PRECISION=DOUBLE required)')
    if df['SINGLE_PRECISION_JACOBIAN']:
        die('SINGLE_PRECISION_JACOBIAN (jac_t = float) is not supported by the kernel')
    if df['NONAKA_PLOT']:
        die('NONAKA_PLOT changes burn_t (level, reference_time, strang_half): not supported')
    if [v for _, v in d['rates']] != list(range(1, cfg['NumRates'] + 1)):
        die(f'Rates enumerators are not 1..NumRates in order: {d["rates"]}')
    if [v for _, v in d['species']] != list(range(1, cfg['NumSpecTotal'] + 1)):
        die(f'Species enumerators are not 1..NumSpecTotal in order: {d["species"]}')
    for s, (a_f, z_f, a_arr, z_arr, _) in d['ion'].items():
        if a_f != a_arr or z_f != z_arr:
            die(f'species {s}: NetworkProperties::aion/zion differ from aion[]/zion[]')
    for r_, (fwd, rev) in d['densexp'].items():
        if fwd > 3 or rev > 3:
            die(f'rate {r_}: density exponent > 3 (rhs.H static_asserts this too)')
    for r_, rdv in d['rd'].items():
        if max(rdv['exponent']) > 3:
            die(f'rate {r_}: species exponent > 3 is not supported by rhs_term')
    if df['SCREENING'] and d['screen_method'] not in SCREEN_FIELDS:
        die(f'SCREEN_METHOD {d["screen_method"]} not supported')
    warn = []
    if df['NEUTRINOS']:
        warn.append('NEUTRINOS is defined: network_rhs.H refuses to compile (sneut5 is not ported); '
                    'rebuild the flags with USE_NEUTRINOS=FALSE for a kernel today')
    if df['NSE'] or df['NSE_NET']:
        warn.append('NSE/NSE_NET is defined: burn.H carries the extra fields, but the NSE burner '
                    'path is not ported (network_rhs.H refuses to compile)')
    if not (df['STRANG']):
        warn.append('not a STRANG build: the ported integrator glue is Strang only')

    gen = Gen(d, cmdline, flags_src)
    texts = {'network_data.H': gen.network_data(), 'burn.H': gen.burn()}

    # verify in a scratch copy, then install
    with tempfile.TemporaryDirectory(dir=out) as tmp:
        for n, t in texts.items():
            open(os.path.join(tmp, n), 'w').write(t)
        chk = os.path.join(out, 'check_generated.cpp')
        open(chk, 'w').write(gen.check_tu(tmp))
        cexe = os.path.join(out, 'check_generated')
        r = run([a.cxx, *flags, f'-I{a.support_dir}', '-O0', '-w', chk, '-o', cexe], cwd=test)
        open(os.path.join(out, 'check_build.log'), 'w').write(r.stdout + r.stderr)
        if r.returncode:
            die(f'generated headers disagree with the stock source, see {out}/check_build.log:\n' +
                '\n'.join(l for l in r.stderr.splitlines() if 'error' in l)[:3000])
        r = run([cexe], cwd=test)
        if r.returncode:
            die(f'{r.returncode} aion/zion/aion_inv array entries differ from the generated lists')
        tu = os.path.join(tmp, 'standalone.cpp')
        open(tu, 'w').write('#include "network_data.H"\n#include "burn.H"\nint main() { return 0; }\n')
        r = run([a.cxx, '-std=c++20', '-fsyntax-only', '-Wall', '-Werror', f'-I{a.support_dir}', tu])
        if r.returncode:
            die('standalone syntax check failed:\n' + r.stderr[:3000])
        for n in texts:
            shutil.copyfile(os.path.join(tmp, n), os.path.join(emit, n))

    pairs = len(gen.screen_pairs()) if df['SCREENING'] else 0
    print(f'gen_network_data: {d["network"]}: NumSpec={cfg["NumSpec"]} NumSpecTotal={cfg["NumSpecTotal"]} '
          f'NumAux={cfg["NumAux"]} NumRates={cfg["NumRates"]} intermediate={cfg["num_intermediate"]} '
          f'screen_pairs={pairs} bion/mion={"yes" if have_bm else "no"}')
    print(f'gen_network_data: wrote {os.path.join(emit, "network_data.H")} and {os.path.join(emit, "burn.H")} '
          '(checked against the stock source)')
    for w in warn:
        print('gen_network_data: WARNING: ' + w)


if __name__ == '__main__':
    main()

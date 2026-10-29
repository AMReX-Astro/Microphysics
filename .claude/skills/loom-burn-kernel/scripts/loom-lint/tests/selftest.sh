#!/bin/bash
# Lint self-test: every rule must parse and fire on should_fire.cpp.
# Run after adding or editing a rule (a rule with a YAML error makes ast-grep
# skip the whole scan with "Error: Cannot parse rule", which is easy to miss).
set -u
here=$(cd "$(dirname "$0")" && pwd)
cfg="$here/../sgconfig.yml"
err=$(mktemp)
out=$(ast-grep scan -c "$cfg" --json=compact "$here/should_fire.cpp" 2> "$err")
if grep -q "Cannot parse rule" "$err"; then echo "FAIL:"; head -8 "$err"; rm -f "$err"; exit 1; fi
rm -f "$err"
fired=$(echo "$out" | python3 -c 'import json,sys; print("\n".join(sorted({d["ruleId"] for d in json.load(sys.stdin)})))')
status=0
for r in "$here"/../rules/*.yml; do
  id=$(basename "$r" .yml)
  if echo "$fired" | grep -qx "$id"; then echo "ok    $id"; else echo "MISS  $id"; status=1; fi
done
[ $status -eq 0 ] && echo PASS || echo FAIL
exit $status

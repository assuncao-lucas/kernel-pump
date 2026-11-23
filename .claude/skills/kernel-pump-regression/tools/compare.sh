#!/usr/bin/env bash
# Compare two capture dirs field by field, ignoring the two wall-clock fields.
# Usage: compare.sh <dir-A> <dir-B>
set -u
A="$1"; B="$2"
strip() { grep -vE '^(total time \(s\)|time building kernel and buckets \(s\)):' "$1"; }
same=0; differ=0; missing=0
for fa in "$A"/*.sol; do
  f=$(basename "$fa"); fb="$B/$f"
  if [ ! -f "$fb" ]; then echo "MISSING in B: $f"; missing=$((missing+1)); continue; fi
  if diff -q <(strip "$fa") <(strip "$fb") >/dev/null; then
    same=$((same+1))
  else
    differ=$((differ+1)); echo "=== DIFFERS: $f"; diff <(strip "$fa") <(strip "$fb") | sed 's/^/    /'
  fi
done
for fb in "$B"/*.sol; do f=$(basename "$fb"); [ -f "$A/$f" ] || { echo "MISSING in A: $f"; missing=$((missing+1)); }; done
echo "----"
echo "identical: $same   differing: $differ   missing: $missing"
[ "$differ" -eq 0 ] && [ "$missing" -eq 0 ]

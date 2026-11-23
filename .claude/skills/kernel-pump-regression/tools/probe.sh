#!/usr/bin/env bash
# Single-run probe: find instances that exercise the KP bucket loop but finish inside the budget.
set -u
REPO=/home/lucas/Documents/Research/kernel-pump
KP="$1"; INSTDIR="$2"; OUT="$3"; TL="$4"; CFG="$5"; shift 5
mkdir -p "$OUT/sol"
printf 'instance\ttime\tstatus\titers\tbuckets\tlastbkt\tfirstiter\n' > "$OUT/probe.tsv"
for inst in "$@"; do
  base="${inst%.mps}"; base="${base%.gz}"
  timeout $((TL + 120)) "$KP" "$INSTDIR/$inst" -c "$REPO/settings/all/$CFG.cfg" \
      seed=1 solutionFolder="$OUT/sol" timeLimit="$TL" >/dev/null 2>&1
  f="$OUT/sol/s_${CFG}_${base}_1.sol"
  if [ -f "$f" ]; then
    t=$(grep -oP '^total time \(s\): \K.*' "$f")
    st=$(grep -oP '^STATUS: \K.*' "$f" | sed 's/FAILED.*/FAIL/;s/FOUND INTEGER FEASIBLE/FOUND/')
    it=$(grep -oP '^# iterations: \K.*' "$f"); nb=$(grep -oP '^# buckets: \K.*' "$f")
    lb=$(grep -oP '^last bucket visited: \K.*' "$f"); fi_=$(grep -oP '^first bucket to iter pump: \K.*' "$f")
  else t=-; st=NOFILE; it=-; nb=-; lb=-; fi_=-; fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$base" "$t" "$st" "$it" "$nb" "$lb" "$fi_" >> "$OUT/probe.tsv"
  echo "[$base] ${t}s $st iters=$it buckets=$nb lastbkt=$lb"
done

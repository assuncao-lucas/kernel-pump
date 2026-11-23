#!/usr/bin/env bash
# M0 determinism screen: run each (instance,config,seed) twice with the SAME time limit
# and compare .sol files ignoring the two wall-clock fields.
set -u
REPO=/home/lucas/Documents/Research/kernel-pump
KP="$1"; INSTDIR="$2"; OUTROOT="$3"; TL="$4"; CONFIGS="$5"; SEEDS="$6"; INSTANCES="$7"

# fields that legitimately vary run to run
strip_timing() { grep -vE '^(total time \(s\)|time building kernel and buckets \(s\)):' "$1"; }

mkdir -p "$OUTROOT/A" "$OUTROOT/B"
printf 'config\tinstance\tseed\tdeterministic\ttimeA\ttimeB\tstatusA\titersA\tlastbktA\n' > "$OUTROOT/screen.tsv"

for cfg in $CONFIGS; do
  for inst in $INSTANCES; do
    base="${inst%.mps}"
    for seed in $SEEDS; do
      for run in A B; do
        "$KP" "$INSTDIR/$inst" -c "$REPO/settings/all/$cfg.cfg" \
          seed="$seed" solutionFolder="$OUTROOT/$run" timeLimit="$TL" >/dev/null 2>&1
      done
      fa="$OUTROOT/A/s_${cfg}_${base}_${seed}.sol"
      fb="$OUTROOT/B/s_${cfg}_${base}_${seed}.sol"
      if [ -f "$fa" ] && [ -f "$fb" ]; then
        if diff -q <(strip_timing "$fa") <(strip_timing "$fb") >/dev/null; then det=YES; else det=NO; fi
        ta=$(grep -oP '^total time \(s\): \K.*' "$fa"); tb=$(grep -oP '^total time \(s\): \K.*' "$fb")
        st=$(grep -oP '^STATUS: \K.*' "$fa" | sed 's/FAILED TO FIND AN INTEGER FEASIBLE SOLUTION/FAIL/;s/FOUND INTEGER FEASIBLE/FOUND/')
        it=$(grep -oP '^# iterations: \K.*' "$fa"); lb=$(grep -oP '^last bucket visited: \K.*' "$fa")
      else
        det=NOFILE; ta=-; tb=-; st=-; it=-; lb=-
      fi
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cfg" "$base" "$seed" "$det" "$ta" "$tb" "$st" "$it" "$lb" >> "$OUTROOT/screen.tsv"
      echo "[$cfg $base seed=$seed] $det  (A=${ta}s B=${tb}s $st)"
    done
  done
done

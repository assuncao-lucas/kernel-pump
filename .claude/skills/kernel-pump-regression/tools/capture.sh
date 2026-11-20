#!/usr/bin/env bash
# Capture .sol files for a (configs x instances x seeds) grid with one binary.
# Usage: capture.sh <kp-binary> <instance-dir> <out-dir> <time-limit> <configs-csv> <seeds-csv> <instances-csv>
set -u
REPO=/home/lucas/Documents/Research/kernel-pump
KP="$1"; INSTDIR="$2"; OUT="$3"; TL="$4"
CONFIGS=$(echo "$5" | tr ',' ' '); SEEDS=$(echo "$6" | tr ',' ' '); INSTANCES=$(echo "$7" | tr ',' ' ')
mkdir -p "$OUT"
{
  echo "binary:      $KP"
  echo "binary_sha:  $(sha256sum "$KP" | cut -c1-16)"
  echo "commit:      $(git -C "$REPO" rev-parse HEAD)"
  echo "dirty:       $(git -C "$REPO" status --porcelain -- src include CMakeLists.txt | tr '\n' ';')"
  echo "instancedir: $INSTDIR"
  echo "timelimit:   $TL"
  echo "configs:     $5"
  echo "seeds:       $6"
  echo "instances:   $7"
  echo "captured:    $(date -Is)"
} > "$OUT/MANIFEST"
n=0
for cfg in $CONFIGS; do for inst in $INSTANCES; do for seed in $SEEDS; do
  n=$((n+1))
  "$KP" "$INSTDIR/$inst.mps" -c "$REPO/settings/all/$cfg.cfg" \
      seed="$seed" solutionFolder="$OUT" timeLimit="$TL" >/dev/null 2>&1
  echo "[$n] $cfg $inst seed=$seed -> $(grep -oP '^STATUS: \K.*' "$OUT/s_${cfg}_${inst}_${seed}.sol" 2>/dev/null | head -c 30)"
done; done; done
echo "captured $n runs into $OUT"

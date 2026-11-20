#!/usr/bin/env bash
#
# Is a case reproducible at all? Run it N times with ONE binary and count how many distinct
# behavioural fingerprints come back. 1 = reproducible. More than 1 = the case is flaky and cannot
# be used as a regression reference, whatever a verify.sh difference suggested.
#
# Usage:
#   stability.sh <binary> <instance-dir> <out-dir> <time-limit> <N> <instance> [configs...]
#
# Example -- is air05 reproducible under the four KP configs?
#   tools/stability.sh golden/bin/kp_baseline /path/to/instances /tmp/stab 30 3 air05 \
#       config16 config17 config21 config23
#
# Run it on an idle machine: the thing being measured is sensitivity to timing.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"

KP="${1:?usage: stability.sh <binary> <instance-dir> <out-dir> <time-limit> <N> <instance> [configs...]}"
INSTDIR="${2:?instance dir}"
OUT="${3:?out dir}"
TL="${4:?time limit}"
N="${5:?number of runs}"
INST="${6:?instance (without .mps)}"
shift 6
CONFIGS=("$@")
[ ${#CONFIGS[@]} -eq 0 ] && CONFIGS=(config16 config17 config21 config23)

mkdir -p "$OUT"
printf 'config\tinstance\tseed\tdistinct_of_%s\tfingerprints\n' "$N" >"$OUT/stability.tsv"

flaky=0
for cfg in "${CONFIGS[@]}"; do
	for seed in 1 2 3; do
		sigs=""
		for _ in $(seq 1 "$N"); do
			rm -rf "$OUT/r"
			mkdir -p "$OUT/r"
			"$KP" "$INSTDIR/$INST.mps" -c "$REPO/settings/all/$cfg.cfg" seed="$seed" \
				solutionFolder="$OUT/r" timeLimit="$TL" >/dev/null 2>&1
			f="$OUT/r/s_${cfg}_${INST}_${seed}.sol"
			if [ -f "$f" ]; then
				s=$(grep -vE '^(total time \(s\)|time building kernel and buckets \(s\)):' "$f" |
					md5sum | cut -c1-8)
			else
				s=NOFILE
			fi
			sigs="$sigs $s"
		done
		d=$(echo $sigs | tr ' ' '\n' | sort -u | grep -c .)
		[ "$d" -gt 1 ] && flaky=$((flaky + 1))
		printf '%s\t%s\t%s\t%s\t[%s]\n' "$cfg" "$INST" "$seed" "$d" "$sigs" >>"$OUT/stability.tsv"
		printf '%-10s seed=%s  distinct=%s :%s\n' "$cfg" "$seed" "$d" "$sigs"
	done
done
rm -rf "$OUT/r"

echo "----"
if [ "$flaky" -eq 0 ]; then
	echo "every case reproduced identically over $N runs"
else
	echo "$flaky case(s) are NOT reproducible over $N runs -- they belong in golden/EXCLUSIONS.tsv,"
	echo "not in a bug report. See golden/stability-neos-1456979.tsv for the precedent."
fi
echo "detail: $OUT/stability.tsv"

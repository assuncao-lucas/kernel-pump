#!/usr/bin/env bash
#
# Decide what a verify.sh difference actually means. For each case it runs BOTH the baseline binary
# and the candidate binary N times and compares the three things that matter: what the candidate
# does, what the baseline does now, and what the golden file recorded.
#
# Verdicts:
#   REGRESSION     both binaries are self-consistent but disagree with each other.
#                  This is the only verdict that means your change broke something.
#   FLAKY          the baseline disagrees with itself over N runs. The case cannot serve as a
#                  reference at all; it belongs in golden/EXCLUSIONS.tsv.
#   STALE GOLDEN   both binaries agree with each other but not with the golden file. The recorded
#                  value is no longer reproducible on this machine -- re-record it from the
#                  baseline, or drop the case.
#   OK             everything agrees; the difference verify.sh saw was transient.
#
# Usage:
#   adjudicate.sh <candidate-binary> <instance-dir> <N> <config:instance:seed:timelimit> ...
#
# Example (the two cases verify.sh flagged):
#   tools/adjudicate.sh ../../../build/kp /path/to/instances 3 \
#       config21:air05:2:30 config21:air05:3:30
#
# Run it on an idle machine.
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL="$(cd "$HERE/.." && pwd)"
REPO="$(cd "$SKILL/../../.." && pwd)"
GOLDEN="$SKILL/golden"
BASELINE="$GOLDEN/bin/kp_baseline"

CAND="${1:?usage: adjudicate.sh <candidate-binary> <instance-dir> <N> <config:instance:seed:tl> ...}"
INSTDIR="${2:?instance dir}"
N="${3:?number of runs}"
shift 3

if [ ! -x "$BASELINE" ]; then
	echo "adjudicate.sh: baseline binary missing at $BASELINE" >&2
	echo "               it is git-ignored; rebuild it from $(cat "$GOLDEN/bin/BASELINE_COMMIT" 2>/dev/null)" >&2
	exit 2
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/kp-adjudicate-XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

strip() { grep -vE '^(total time \(s\)|time building kernel and buckets \(s\)):' "$1"; }

# fingerprints of N runs of one binary, space separated
fingerprints() {
	local bin="$1" cfg="$2" inst="$3" seed="$4" tl="$5" out
	local sigs=""
	for _ in $(seq 1 "$N"); do
		out="$WORK/run"
		rm -rf "$out"
		mkdir -p "$out"
		"$bin" "$INSTDIR/$inst.mps" -c "$REPO/settings/all/$cfg.cfg" seed="$seed" \
			solutionFolder="$out" timeLimit="$tl" >/dev/null 2>&1
		local f="$out/s_${cfg}_${inst}_${seed}.sol"
		if [ -f "$f" ]; then
			sigs="$sigs $(strip "$f" | md5sum | cut -c1-8)"
		else
			sigs="$sigs NOFILE"
		fi
	done
	echo "$sigs"
}

uniq_of() { echo "$1" | tr ' ' '\n' | grep -c . >/dev/null; echo "$1" | tr ' ' '\n' | sed '/^$/d' | sort -u; }

regressions=0
for spec in "$@"; do
	IFS=':' read -r cfg inst seed tl <<<"$spec"
	: "${cfg:?}" "${inst:?}" "${seed:?}" "${tl:?}"

	gf="$GOLDEN/sol/s_${cfg}_${inst}_${seed}.sol"
	gold=$([ -f "$gf" ] && strip "$gf" | md5sum | cut -c1-8 || echo MISSING)

	base_sigs=$(fingerprints "$BASELINE" "$cfg" "$inst" "$seed" "$tl")
	cand_sigs=$(fingerprints "$CAND" "$cfg" "$inst" "$seed" "$tl")
	base_set=$(uniq_of "$base_sigs" | tr '\n' ',' | sed 's/,$//')
	cand_set=$(uniq_of "$cand_sigs" | tr '\n' ',' | sed 's/,$//')

	echo "=== $cfg / $inst / seed=$seed / tl=$tl"
	echo "    golden   : $gold"
	echo "    baseline :$base_sigs"
	echo "    candidate:$cand_sigs"

	if [ "$(uniq_of "$base_sigs" | grep -c .)" -gt 1 ]; then
		echo "    VERDICT  : FLAKY -- the baseline disagrees with itself; exclude this case"
	elif [ "$base_set" != "$cand_set" ]; then
		echo "    VERDICT  : REGRESSION -- both binaries are self-consistent but disagree"
		regressions=$((regressions + 1))
	elif [ "$base_set" != "$gold" ]; then
		echo "    VERDICT  : STALE GOLDEN -- both binaries agree ($base_set) but the record says $gold;"
		echo "               re-record from the baseline, or drop the case"
	else
		echo "    VERDICT  : OK -- everything agrees; the verify.sh difference was transient"
	fi
	echo
done

echo "----"
if [ "$regressions" -eq 0 ]; then
	echo "no regression found: nothing here is attributable to the candidate binary"
else
	echo "$regressions case(s) look like a REAL regression in the candidate"
fi
[ "$regressions" -eq 0 ]

#!/usr/bin/env bash
#
# Replay the golden regression cases with a candidate kp binary and compare, field by field,
# against the recorded baseline output. Exits non-zero if ANY case differs behaviourally.
#
# Usage:
#   ./verify.sh                          # tests <repo>/build/kp
#   ./verify.sh /path/to/kp              # tests a specific binary
#   ./verify.sh /path/to/kp /path/to/instances
#
# The instance directory may also be given as $KP_INSTANCE_DIR. It must contain the MIPLIB .mps
# files named in golden/CASES.tsv (the same files instances/miplib.txt lists).
#
# Ignored when comparing: 'total time (s)' and 'time building kernel and buckets (s)'. Everything
# else must match exactly -- STATUS, # iterations, # buckets, last bucket visited, first bucket to
# iter pump, value, reopt value, both integrality gaps, num frac, and all four 'num bin vars' lines.
#
# Run it on an idle machine. KP is not unconditionally reproducible: the per-bucket budget is
# time_left/(num_buckets+1), derived from the wall clock, so CPU contention can legitimately change
# a bucket's iteration count. Never compile while this is running.

set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GOLDEN="$HERE/golden"
REPO="$(cd "$HERE/../../.." && pwd)"

KP="${1:-$REPO/build/kp}"
INSTDIR="${2:-${KP_INSTANCE_DIR:-/home/lucas/Downloads/instances/benchmark}}"

if [ ! -x "$KP" ]; then
	echo "verify.sh: no kp binary at $KP" >&2
	echo "           build one first, or pass the path as the first argument" >&2
	exit 2
fi
if [ ! -d "$INSTDIR" ]; then
	echo "verify.sh: instance directory not found: $INSTDIR" >&2
	echo "           pass it as the second argument or set KP_INSTANCE_DIR" >&2
	exit 2
fi
if [ ! -f "$GOLDEN/CASES.tsv" ]; then
	echo "verify.sh: $GOLDEN/CASES.tsv is missing" >&2
	exit 2
fi

# Scratch output goes outside the repo so a verification run can never dirty the working tree.
NEW="$(mktemp -d "${TMPDIR:-/tmp}/kp-verify-XXXXXX")"
trap 'rm -rf "$NEW"' EXIT

FAILLOG="$NEW/VERIFY_FAILURES.txt"
: >"$FAILLOG"

strip() { grep -vE '^(total time \(s\)|time building kernel and buckets \(s\)):' "$1"; }

echo "candidate : $KP"
echo "            sha256(16) $(sha256sum "$KP" | cut -c1-16)"
echo "golden    : $GOLDEN  ($(($(wc -l <"$GOLDEN/CASES.tsv") - 1)) cases)"
echo "instances : $INSTDIR"
echo

same=0
differ=0
n=0
total=$(($(wc -l <"$GOLDEN/CASES.tsv") - 1))

while IFS=$'\t' read -r tier cfg inst seed tl; do
	[ "$tier" = "tier" ] && continue
	n=$((n + 1))
	"$KP" "$INSTDIR/$inst.mps" -c "$REPO/settings/all/$cfg.cfg" \
		seed="$seed" solutionFolder="$NEW" timeLimit="$tl" >/dev/null 2>&1
	f="s_${cfg}_${inst}_${seed}.sol"
	if [ ! -f "$NEW/$f" ]; then
		differ=$((differ + 1))
		echo "NO OUTPUT: $f (tier=$tier tl=$tl)" | tee -a "$FAILLOG"
	elif diff -q <(strip "$GOLDEN/sol/$f") <(strip "$NEW/$f") >/dev/null; then
		same=$((same + 1))
	else
		differ=$((differ + 1))
		{
			echo "=== DIFFERS: $f  (tier=$tier tl=$tl)"
			diff <(strip "$GOLDEN/sol/$f") <(strip "$NEW/$f") | sed 's/^/    /'
		} | tee -a "$FAILLOG"
	fi
	printf '\r[%d/%d] identical=%d differing=%d' "$n" "$total" "$same" "$differ"
done <"$GOLDEN/CASES.tsv"

echo
echo "----"
echo "golden cases: $total   identical: $same   differing: $differ"
if [ "$differ" -eq 0 ]; then
	echo "VERDICT: behaviour unchanged"
else
	echo "VERDICT: REGRESSION in $differ case(s)"
	# kept outside the repo so a verification run never dirties the working tree
	keep="${TMPDIR:-/tmp}/kp-verify-failures-$(date +%Y%m%d-%H%M%S).txt"
	cp "$FAILLOG" "$keep" 2>/dev/null && echo "detail: $keep"
	echo
	echo "A difference is NOT necessarily a regression. Adjudicate it -- this runs the baseline and"
	echo "the candidate 3 times each and tells you whether it is a REGRESSION, a FLAKY case or a"
	echo "STALE GOLDEN record:"
	echo
	echo "  $HERE/tools/adjudicate.sh $KP $INSTDIR 3 <config>:<instance>:<seed>:<timelimit>"
	echo
	echo "e.g. config21:air05:2:30 . See the skill's SKILL.md for what each verdict means."
fi

[ "$differ" -eq 0 ]

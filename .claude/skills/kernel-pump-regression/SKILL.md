---
name: kernel-pump-regression
description: Behavioural regression suite for this Feasibility Pump / Kernel Pump codebase — a golden set of 191 recorded .sol outputs plus verify.sh, which replays them with a candidate kp binary and reports any behavioural difference. Use after changing anything that can affect FP or KP results (src/, include/kernelpump/, extern/, settings/), before committing or merging such a change, when asked whether a refactor preserved behaviour, or when a run's status, iteration count, bucket count or integrality gap looks unexpected. Also explains why KP is not bit-reproducible and how to tell a real regression from a flaky case.
---

# Kernel Pump behavioural regression suite

The project has no general unit-test suite, and FP/KP correctness is mostly emergent: a refactor
that compiles, runs and produces a feasible solution can still have changed the search trajectory
and therefore every number in the article. This suite is the safety net for that.

It answers one question: **does a candidate `kp` binary still behave exactly like the recorded
baseline?**

## Running it

```bash
.claude/skills/kernel-pump-regression/verify.sh                   # tests ./build/kp
.claude/skills/kernel-pump-regression/verify.sh /path/to/kp       # tests a specific binary
.claude/skills/kernel-pump-regression/verify.sh /path/to/kp /path/to/instances
```

The instance directory defaults to `/home/lucas/Downloads/instances/benchmark` and can also be set
with `$KP_INSTANCE_DIR`. It must hold the MIPLIB `.mps` files named in `golden/CASES.tsv`.

Exit code 0 means every case matched. Non-zero names the differing cases and writes the detail to a
file outside the repo. Expect **roughly 30 minutes**; it replays 191 CPLEX runs sequentially.

Two rules that are not optional:

- **Run it on an idle machine, and never compile while it runs.** The per-bucket time budget is
  `time_left / (num_buckets + 1)`, taken from the wall clock, so CPU contention can legitimately
  change a bucket's iteration count. A parallel build once produced two phantom "regressions" that
  vanished on a quiet machine.
- **Test a fixed binary.** Rebuilding `build/kp` underneath a running verification compares against
  a moving target.

## What counts as a difference

Every field of the `.sol` record except the two wall-clock ones (`total time (s)` and
`time building kernel and buckets (s)`). So: `STATUS`, `# iterations`, `# buckets`,
`last bucket visited`, `first bucket to iter pump`, `value`, `reopt value`,
`real integrality gap`, `projection integrality gap`, `num frac`, and all four `num bin vars *`
lines.

The two sharpest detectors, worth reading first in any failure:

- **`# iterations`** catches an accidental reset of FP's shared iteration budget. `fp.iterLimit` is
  decremented as iterations are consumed and is *not* restored by `init()`, so it is a budget
  shared by every bucket. Anything that gives each bucket a fresh budget — most obviously
  constructing a new heuristic per bucket — shows up here immediately.
- **`last bucket visited`** catches any change in the kernel-absorption branches.

## Contents

| Path | What it is |
| --- | --- |
| `verify.sh` | the entry point described above |
| `golden/MANIFEST` | how the set was captured: baseline commit, compiler, CPLEX version, tiers, the determinism screen and its limitations |
| `golden/CASES.tsv` | the authoritative case list — tier, config, instance, seed, time limit |
| `golden/sol/` | 191 recorded baseline `.sol` files, one per case |
| `golden/EXCLUSIONS.tsv` | cases deliberately left out, each with its reason |
| `golden/stability-neos-1456979.tsv` | the evidence that dropped one instance (see below) |
| `golden/bin/BASELINE_COMMIT` | the commit the golden output came from |
| `golden/bin/kp_baseline` | the baseline binary — **git-ignored on purpose** (34 MB). Rebuild it from `BASELINE_COMMIT` if you need it |
| `tools/screen.sh` | run a grid twice and emit a per-case `deterministic YES/NO` verdict |
| `tools/probe.sh` | single-run probe reporting time, status, iterations, buckets — used to pick instances |
| `tools/capture.sh` | capture a grid into a directory with a `MANIFEST` |
| `tools/compare.sh` | diff two capture directories field by field |

## How the set is built, and why it looks like this

Three tiers, chosen so the expensive instances are not needlessly repeated across seeds:

| tier | content | grid | time limit |
| --- | --- | --- | --- |
| A | fast instances, the workhorse | 4 KP configs × 13 instances × 3 seeds, minus 1 flaky case = **155** | 30s |
| B | slow instances with maximal bucket coverage | 4 KP configs × 4 instances × 1 seed = **16** | 120s |
| C | the FP and plain-solver paths | 5 configs × 4 instances × 1 seed = **20** | 30s |

KP configs are `config16`, `config17`, `config21`, `config23`; FP are `config1`, `config3`,
`config20`; plain solver are `config18`, `config19` — the nine of `settings/settings_list.txt`.
Tier C exists so a change aimed at KP cannot quietly disturb `method=feaspump` or `method=solver`.

**Do not pick instances by size.** Of the 20 smallest instances in `instances/miplib.txt`, 17 report
`last bucket visited: 0` with 0–1 iterations: KP solves them on the initial kernel and never enters
the bucket loop, so they cannot detect a regression in the interesting code. Select on behaviour
instead — `last bucket visited > 0`, and a total time comfortably under the limit. `tools/probe.sh`
exists for that. A run that exhausts `fp.iterLimit` is ideal: it terminates early *and*
deterministically, since once the shared budget is gone the remaining buckets get zero iterations.
94 of the 191 cases exercise the bucket loop; `traininstance2` walks all 5 of its buckets
with 11402 iterations, `ns1830653` all 11, `csched007` all 8.

## KP is not bit-reproducible — telling flaky from broken

Because the per-bucket budget comes from the wall clock, two runs of the same
instance × config × seed can legitimately differ even when neither hits its time limit. A case that
finishes in 5s under a 30s limit is *not* automatically safe.

So when `verify.sh` reports a difference, **do not assume a regression**. First establish whether
the case is even self-consistent on the baseline:

```bash
.claude/skills/kernel-pump-regression/tools/screen.sh \
    .claude/skills/kernel-pump-regression/golden/bin/kp_baseline \
    /path/to/instances /tmp/screen 30 "config17" "1" "someinstance.mps"
```

If the baseline disagrees with itself, the case is flaky and belongs in `EXCLUSIONS.tsv`, not in a
bug report.

`golden/stability-neos-1456979.tsv` is the cautionary example. That instance **passed** a two-run
determinism screen on several cases, yet a three-run repeat of the *unmodified baseline* produced up
to **3 distinct outcomes per case** (config17: 3/3 distinct for every seed). All 8 of its cases were
removed, taking the set from 199 to 191.

**If you rebuild this set, screen with at least 3 runs, not 2.** Two agreeing runs are not a
determinism proof.

## Related

- `docs/spec-heuristic-interface.md` — the design this suite was built to validate (the pluggable
  `FeasibilityHeuristic` interface), its lifecycle contract and its behaviour-preservation traps.
- `tests/heuristic_interface_test.cpp` — the complementary *mechanism* test: golden files prove the
  outcome, that proves KP calls the heuristic interface in the promised order and covers the
  kernel-absorption branches no config in `settings/all` reaches.
- The `kernel-pump-internals` skill — what the algorithms actually do, and how to read a `.sol`.

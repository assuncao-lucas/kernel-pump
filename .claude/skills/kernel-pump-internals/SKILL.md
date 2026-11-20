---
name: kernel-pump-internals
description: Full behavioural reference for this Feasibility Pump / Kernel Pump MIP heuristic codebase — the algorithms, the solver-abstraction layer, how kernel/bucket construction works, how the FP and KP control flow and budgets interact, config key semantics, and the non-obvious invariants that break things silently. Use when reading, modifying, debugging, extending or reviewing any C++ code under src/, include/kernelpump/, extern/, or results/, when interpreting .sol output or summary.csv, when changing config files in settings/, or when a run produces an unexpected status, iteration count, or integrality gap.
---

# Kernel Pump / Feasibility Pump internals

A C++17 research codebase implementing the **Feasibility Pump (FP)** family of MIP primal
heuristics plus a novel variant, the **Kernel Pump (KP)**. Both look for a *feasible* solution to
a MIP — neither is an optimiser. Built by extending the public FeasPumpCollection.

Everything below was derived by reading the implementation. Where behaviour is surprising, the
file:line is given so you can re-verify rather than trust this document.

## The one thing to understand first

A "kernel" is **not** a submodel. Variables are never added or removed. KP works entirely by
moving the **upper bound of binary variables between 0 and 1**:

- *activate* a binary variable → `ub := 1` (free to take either value)
- *deactivate* → `ub := 0` (pinned to 0)

See `CPXModel::updateModelVarBounds` (`src/cpxmodel.cpp:972`). The column count never changes.
So "the kernel" = the set of binaries currently allowed to be nonzero; every other binary is
clamped to 0. Any reasoning that assumes a shrinking model is wrong.

## Reading map

| Concern | Where |
| --- | --- |
| CLI entry, method dispatch, final feasibility re-check, output | `src/main.cpp` |
| FP algorithm | `src/feaspump.cpp` (2.5k lines), `include/kernelpump/feaspump.h` |
| KP algorithm | `src/kernelpump.cpp`, `include/kernelpump/kernelpump.h` |
| Solver abstraction (and KP-specific extensions) | `include/kernelpump/mipmodel.h` |
| Only working backend | `src/cpxmodel.cpp` (CPLEX) |
| Rounding strategies ("transformers") | `src/transformers.cpp` |
| Variable orderings for rounding ("rankers") | `src/ranking.cpp` |
| Output record | `src/solution.cpp`, `include/kernelpump/solution.h` |
| Results aggregation (separate binary) | `results/src/main.cpp` |
| Config key semantics, config№ → paper name | `references/config-reference.md` |

`src/xprsmodel.cpp`, `src/scipmodel.cpp`, `src/pdlpmodel.cpp` exist and compile conditionally but
are **not production-ready**. Only `solver=cpx` works. Don't spend effort on them unless asked.

## Three methods, one entry point

`method=` in the config selects the path in `src/main.cpp` (case-insensitive; empty or unrecognised
throws, and the error message lists the registered heuristics):

- **`solver`** — plain CPLEX MIP solve with `SolutionLimit=1`, i.e. stop at the first feasible
  solution. The baseline. `mipFeasEmphasis=1` additionally sets CPLEX feasibility emphasis.
- **`kernelpump`** — `KernelPump`, which drives a feasibility heuristic repeatedly over a growing
  kernel (which heuristic: `kp.heuristic`, default `feaspump`).
- **any registered feasibility-heuristic name** — run that heuristic standalone over the whole
  problem. `feaspump` is the only one registered today, so in practice `method=feaspump` means one
  `FeasibilityPump` over the whole problem, exactly as before.

`main.cpp` holds no reference to any concrete heuristic: it resolves the name through
`resolveHeuristicName()` against `HeuristicsFactory`. Registering a new heuristic therefore makes it
available as both `method=<name>` and `kp.heuristic=<name>` with no change to `main.cpp` or
`KernelPump`. See **Pluggable feasibility heuristic** below.

After any of them finds `x`, `main.cpp` does three things that matter:
1. **Re-verifies every constraint** against `x` with tolerance `0.001` and throws if violated
   (`src/main.cpp:334-342`). This is the project's only correctness test — there is no unit suite.
2. **Fixes all integer/binary vars** to their found values, switches to LP, re-solves
   (`src/main.cpp:352-373`). This yields `reopt_value_`: the objective after optimising the
   continuous variables given the integer assignment. `value_` is the raw objective of `x`.
3. Writes `<solutionFolder>/s_<runName>_<probName>_<seed>.sol`.

## FeasibilityPump behaviour

### The loop

`pumpLoop` (`src/feaspump.cpp:1168`) alternates:

1. **round**: `frac_x` → `integer_x` via the configured transformer (`frac2int->apply`).
2. **project**: set the LP objective to the distance to `integer_x`, re-solve the LP → new `frac_x`.
3. stop when `frac_x` is both LP-feasible and integral on the current integer subset.

The objective actually minimised is a blend: `alpha * (original objective, normalised)` +
`(1-alpha) * (distance)`. `runningAlpha` decays by `alphaFactor` each iteration and snaps to 0
below `1e-4` (`src/feaspump.cpp:1963-1965`). With `fp.alpha=0.0` (what the paper's configs use)
it is pure distance minimisation from the start.

### Stages

- **Stage 1** considers **binaries only** (`ignoreGenerals=true`); general integers are ignored in
  both rounding and the distance.
- **Stage 2** considers **all integers**. Entered when stage 1 failed *or* general integers exist
  (`src/feaspump.cpp:765`). It restarts from the best point found in stage 1 (`closestFrac`).
- **Stage 3** (`fp.doStage3`) is a final sub-MIP/RENS attempt. **Off in all paper configs.**

Consequence: on a pure-binary instance that stage 1 solves, stage 2 is skipped entirely.

### Cycle handling

Three escalating mechanisms, all in `pumpLoop`:
- same `integer_x` as last iteration → `perturbe` (flip a few variables).
- `integer_x` already in the cache (`isInCache`) → up to 10 `restart` attempts (stronger random flips).
- counters `pertCnt` / `restartCnt` / `firstPerturbation` are reported in the log.

### Closest point tracking

`closestDist` / `closestPoint` / `closestFrac` record the best (smallest distance) pair seen
(`src/feaspump.cpp:1914-1923`). This is what KP harvests, and what `main.cpp` reports as
`projection_integrality_gap_` when no solution is found.

### Stopping criteria

- hard iteration caps: `fp.iterLimit`, `fp.stage1IterLimit`, `fp.stage2IterLimit`
- time: the `time_limit` argument to `pump()`
- **no-improvement**: if `stopWithNoImprLimit` is on, abort the stage after
  `stage{1,2}NoImprIterLimit` iterations without the distance improving to below
  `imprThreshold ×` the previous best (`src/feaspump.cpp:1914-1972`). With `imprThreshold=0.9`
  that means "40 iterations without a 10% improvement".
- `model->aborted()` (Ctrl-C).

## KernelPump behaviour

### Phase 1 — build kernel and buckets

`BuildKernelAndBuckets` (`src/kernelpump.cpp:183`):

1. Clone the model; optionally **replace the objective** — null objective
   (`buildKernelBasedOnNullObjective`), or `sum of binaries` with MIN or MAX sense
   (`buildKernelBasedOnSumVarsObjective` + `...MaxSense`); optionally flip the sense
   (`reverseObjectiveFunction`).
2. Relax to LP, solve with **dual simplex** (hardcoded `'D'`, chosen so a primal solution is
   available — barrier often yields an optimal value without one, `src/kernelpump.cpp:263`).
   Any of {aborted, failed, infeasible/time-out, not primal feasible} → KP gives up.
3. Sort the binaries. Primary key is the LP **value**, or its **fractional part**
   `|round(v) − v|` when `sortByFractionalPart=1`. Ties break on **reduced cost**. Both keys
   invert for MAX-sense problems (`src/kernelpump.cpp:331-345`).
4. Split into an initial kernel plus an ordered list of buckets, in one of two modes:

**Flat mode** (`buildBucketsByRelaxationLayers=0`): kernel = the first `maxBucketSize` variables;
buckets = consecutive chunks of `maxBucketSize`.

**Layered mode** (`buildBucketsByRelaxationLayers=1`, used by the paper's KP configs): divide the
observed value range into `numBucketLayers` equal-width layers; **re-sort each layer by reduced
cost**; then split each layer further wherever the **sign of the reduced cost** changes. Each
resulting group is a bucket (the first becomes the kernel). Empty groups are dropped, so the
realised bucket count is data-dependent and usually ≠ `numBucketLayers`. Degenerate case: if every
value is identical, group by reduced-cost sign alone (`src/kernelpump.cpp:418-454`).

**Dependency expansion** (`buildBucketsConsideringVariableDependency=1`): when a variable joins a
bucket, also pull in its "dependent" variables — binaries that *co-occur in at least one
constraint* and have LP value > 0 (`addVarToBucket`, `src/kernelpump.cpp:771`). The dependency
relation is the symmetric co-occurrence graph built in `retrieveDependency_impl`
(`include/kernelpump/mipmodel.h:298`) — every pair of variables appearing together in a row. This
is O(nnz²) per row; expensive on dense rows.

**Feasibility repair of the initial kernel** (`tryEnforceFeasibilityInitialKernel=1`): deactivate
all binaries, activate only the kernel, solve the LP. While infeasible, ask CPLEX which
*deactivated* binaries would have to be freed to restore feasibility — via `CPXfeasopt` with
`MIN_SUM`, allowed to relax only the upper bounds of non-kernel binaries (`optimize_set=true`
path of `findSetOfConflictingVariables`, `src/cpxmodel.cpp:1002`) — and absorb those into the
kernel. Repeat until feasible, time-out, or no progress. `num_binary_vars_initial_kernel_` records
the kernel size before this repair, `..._refined_` after.

> An LP-feasible kernel can still be **MIP**-infeasible: the later per-bucket presolve works on
> the integer problem. The code says so explicitly at `src/kernelpump.cpp:682-684`.

### Phase 2 — grow the kernel

`Run` (`src/kernelpump.cpp:798`). All binaries start deactivated. Then iterate
`curr_bucket_index` from `-1` (kernel alone) through the last bucket:

- **Time budget**: `time_left / (num_buckets + 1)` per iteration, except the **last** bucket,
  which gets all remaining time (`src/kernelpump.cpp:826,836`).
- *reference kernel* = current kernel ∪ current bucket. Activate the entering vars, deactivate the
  ones leaving.
- Run `heuristic_->init(model_)` then `->run(...)`. If `init` returns false the subproblem is
  already infeasible and the heuristic never iterates. `heuristic_` is a `FeasibilityHeuristic`
  chosen by `kp.heuristic` (default `feaspump` = the classical FP behind a thin adapter), so KP is
  not tied to FP; see the **Pluggable feasibility heuristic** section below.
- `stopWithNoImprLimit` is **disabled for the last bucket only**, letting the final attempt run
  to the full iteration/time budget (`src/kernelpump.cpp:895`).
- **Warm start**: unless `resetFPBasisAtNewPump=1`, the best fractional point found so far
  (`closest_frac_`) is handed to `pump()` as its starting point, skipping its initial LP.

What happens to the kernel after each FP call — the part worth internalising:

| FP outcome | Effect on the kernel |
| --- | --- |
| **found integer-feasible** | done; absorb reference kernel, record solution, `break` |
| **FP infeasible** (`!feasible_fp`) | absorb the *whole* reference kernel permanently — the bucket's vars stay in forever |
| ran, no solution, `forceBucketVarsIntoKernel=1` | absorb the whole reference kernel |
| ran, no solution, found a strictly better closest point | absorb every binary with `closest_frac > 0`; the rest leave again |
| ran, no solution, no improvement | absorb nothing; all bucket vars leave again |

With `forceBucketVarsIntoKernel=0` the kernel would grow only along variables that the best
fractional point actually uses, and unused bucket variables get clamped back to 0. **But every
config used in the article sets `forceBucketVarsIntoKernel=1`**, so in the published
configuration the kernel simply absorbs each bucket in turn and grows monotonically. The
selective-absorption branches are implemented but unexercised by the paper's settings — verify
against the config in hand before assuming either behaviour.

`last_bucket_visited_` is reported with **0 meaning the initial kernel**.
`first_bucket_to_iter_pump_` is the first bucket where FP managed to iterate at all (i.e. was
LP-feasible and not presolve-infeasible) — useful for diagnosing "FP never got going".

## Budget carry-over — the biggest trap

`KernelPump::Run` calls `heuristic_->readConfig()` **once**, before the loop, but calls `init()` +
`run()` **per bucket**, on a single heuristic instance created in `Init()`. And at the end of every
`pumpLoop`, the iteration limits are **decremented by the iterations just consumed**
(`src/feaspump.cpp:1975-1985`):

```
stage1IterLimit -= used;  stage2IterLimit -= used;  iterLimit -= used;
```

`init()` → `resetPartial()` (`src/feaspump.cpp:450`) does **not** restore them, and does not reset
`nitr`. Therefore:

- **`fp.iterLimit` is a global budget shared across all buckets**, not per bucket. Once exhausted,
  later buckets get zero iterations and fail instantly. A `.sol` reporting exactly
  `# iterations: 12000` against `fp.iterLimit=12000` means the budget ran out, *not* that the last
  bucket worked hard.
- `nitr` accumulates across buckets. Only `resetTotal()` clears it, and KP never calls it.
- Conversely `closestDist`/`closestFrac`/`closestPoint` **are** cleared per bucket by
  `resetPartial`, and `pump()` explicitly re-sets `closestDist = INFBOUND` when given a warm start
  (`src/feaspump.cpp:691`). That is exactly why KP maintains its own `closest_dist_` /
  `closest_frac_` and feeds them back in.

If you change iteration-limit handling, check all three of these interactions.

## Pluggable feasibility heuristic

KP does not depend on FP as such: it needs *some* procedure that, given the sub-MIP induced by the
active kernel, either returns an integer-feasible point or reports the fractional point it got
closest with. That contract is the abstract `FeasibilityHeuristic`
(`include/kernelpump/feasibility_heuristic.h`):

| Member | Role |
| --- | --- |
| `readConfig()` | read own `gConfig()` keys; called **once**, before the first `init()` |
| `init(model, ctype)` | prepare for a subproblem; `false` ⇒ already infeasible, `run()` is skipped |
| `run(timeLimit, HeuristicStart, stopOnNoImprovement)` | the work; returns `HeuristicOutcome{foundSolution, progressed}` |
| `getSolution` / `getPrimalBound` | results on success |
| `getClosestDist` / `getClosestFrac` | best fractional point of the **last** invocation only |
| `getIterations` / `name` | reporting |

`KernelPump` holds one `FeasibilityHeuristicPtr`, created at the end of a successful `Init()` from
`HeuristicsFactory` by the `kp.heuristic` name. `FeasibilityPumpHeuristic`
(`include/kernelpump/feaspump_heuristic.h`, registered in `src/feaspump_heuristic.cpp` as
`"feaspump"`) forwards straight to `FeasibilityPump`, which was left untouched by the
generalisation. Adding a heuristic means implementing the interface and registering it; `KernelPump`
does not change.

Three things that bite:

- **One instance per run.** Implementations may carry state across buckets and FP does — see
  *Budget carry-over* above. Creating one per bucket would hand every bucket a full iteration
  budget and silently invalidate published results.
- **`HeuristicStart::dist` is informational.** `pump()` logs it and then resets its own
  `closestDist` to `INFBOUND` (`src/feaspump.cpp:691`), so a warm start does not seed FP's notion
  of the best distance. KP keeps the running best itself.
- **`progressed` is nearly always true once a warm start exists**, because a warm start assigns
  `closestFrac` immediately (`src/feaspump.cpp:686`) and `progressed` is `!closestFrac.empty()`
  (`src/feaspump.cpp:822`). The `false` branch shows up mainly on the first bucket or with
  `kp.resetFPBasisAtNewPump=1`.

Design, lifecycle contract and the behaviour-preservation traps: `docs/spec-heuristic-interface.md`.
Contract tests: `tests/heuristic_interface_test.cpp` (needs CPLEX), which also covers the
kernel-absorption branches no paper config reaches.

## Presolve discipline

Both FP and KP presolve into a *separate* model and keep the original
(`originalModel` / `original_model_`). `hasPresolve` / `has_presolve_` is set only when presolve
actually reduced something **and** the presolved model still has rows; a presolved model with zero
rows is discarded and presolve undone (`src/feaspump.cpp:530-535`).

Any vector crossing that boundary must be mapped:
- out: `postsolveSolution(x)` (presolved → original)
- in: `presolveSolution(x)` (original → presolved) — KP's warm-start point goes through this on
  every bucket, since each `init()` may presolve afresh (`src/feaspump.cpp:674-682`).

Getting this wrong produces solutions of the wrong length or silently wrong variable values. Always
check which space a vector lives in before using it.

## Rounding and ranking (pluggable, factory-registered)

Transformers (`fp.frac2int`), registered in `src/transformers.cpp:289`:
- **`std`** — `SimpleRounding`: independent threshold rounding; threshold randomised when
  `fp.randomizedRounding=1`.
- **`propround`** — `PropagatorRounding`: rounds variables **one at a time in ranker order**, and
  after each one runs constraint propagation, fixing whatever follows
  (`src/transformers.cpp:217-268`). Far stronger; what every paper config uses. Note
  `fp.filterConstraints` drops numerically nasty rows from the propagator
  (all-continuous with coefficient dynamism > 10, or any row with dynamism > 1000).

Rankers (`fp.ranker`), registered in `src/ranking.cpp:190`:
- **`FRAC`** (default) — most-nearly-integral first, binaries always before general integers
  (score offset `-10`), with `fp.rankNoise` random swaps after `fp.noiseAfter` calls.
- **`LR`** — left to right. **`RND`** — random permutation.

Both the rounder and the ranker are seeded from the global `seed`, which is why `seed` changes
results even though nothing else is stochastic.

## What is live vs. vestigial

`FeasibilityPump` carries a large number of experimental options from prior papers —
multiple reference vectors (`numIntegersObj`, `numFracsObj`, harmonic/exponential-decay weights),
analytic-centre FP (`analcenterFP`), exponential/logistic distance objectives (`expObj`,
`logisObj`), WalkSAT perturbation, and four stage-3 variants (`rensStage3`, `multirensStage3`,
`newStage3`, `normalMIPStage3`). **None are enabled by the configs used in the article.** Before
"fixing" code in those branches, confirm it is actually reachable from the config you care about.
`src/feaspump.cpp` is mostly these branches; the live path is the plain distance objective with a
single reference point.

Also dead: `AddFilesFromDirectory` in `results/src/main.cpp` — every call site is commented out.

## Output

`Solution::WriteToFile` (`src/solution.cpp`) emits one plain-text record. Fields worth knowing:

- `STATUS:` — the only success/failure signal.
- `value` vs `reopt value` — raw objective vs objective after re-optimising continuous vars.
- `real integrality gap` — Σ over integer vars of distance to the nearest value in-domain,
  computed by `computeIntegralityGap` (`include/kernelpump/mipmodel.h:120`) on the **solution** if
  one was found, otherwise on the **closest fractional point**. `num frac` is the count of
  violating variables. On success `num frac` must be 0 — asserted in `main.cpp:253,303`.
- `projection integrality gap` — the FP/KP internal `closestDist`; measured on the possibly
  *presolved* problem and after projection rounding, so it is **not** directly comparable to
  `real integrality gap`. The header says so.
- KP-only fields: `# buckets`, `last bucket visited` (0 = initial kernel),
  `first bucket to iter pump`, `num bin vars added` (final kernel size),
  `num bin vars initial kernel` / `... refined` (before/after feasibility repair).

Batch runs add `summary.csv` (one row per instance × config × seed; `status` and `total time` are
the key columns) and `sorted_files_list.txt` (the instance list actually used — needed to reproduce
a random sample).

## Gotchas

- **`results` rewrites tracked files in place.** `GenerateBoundImprovementLatexTable` overwrites
  `solutions/bound_improvement.txt`. Running the binary to "test" it will dirty the repo's
  committed research output. Check `git status` afterwards.
- **`results` has hardcoded relative paths** (`../../instances/...`, `../../solutions/...`), so it
  only works with CWD = `results/build`. Several of its generator calls are commented out; the
  active one is chosen by editing `ParseArgumentsAndRun`, not by a flag.
- **`-DSILENT_EXEC` is on by default** (top of `CMakeLists.txt`), suppressing most logging. A run
  that "prints nothing" is usually just silent mode. Comment that line out and rebuild to debug.
- **An LP input with no integer variables** makes FP/KP exit immediately having done nothing — they
  have nothing to round. Don't mistake this for a crash. (`build/testando.lp` is such a file.)
- `version.h` is generated at **configure** time, so re-run `cmake`, not just `make`, if you need
  the embedded git hash current.
- Headers must include `<cstdint>`/`<string>` explicitly; newer libstdc++ no longer provides them
  transitively.
- `extern/utils/src/path.cpp` hand-rolls path handling and assumes `/` separators; it is not
  `std::filesystem`.

## Config files

`settings/all/config*.cfg`, listed in `settings/settings_list.txt`. Any key is overridable inline:
`./build/kp inst.mps -c config16.cfg kp.maxBucketSize=20`. For the full key list, the semantics of
each `kp.*` / `fp.*` key, and the config№ → paper-name mapping (FP*, KP*, …), see
`references/config-reference.md`.

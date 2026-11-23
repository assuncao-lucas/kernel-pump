# Spec: Feasibility Jump (`fj`) as a feasibility heuristic, standalone and inside the Kernel Pump

**Status:** specified, not implemented. The open questions were settled on 2026-11-23 (§6.4).
**Branch:** `feasibility-jump`, created from `origin/generic_interface` @ `4c881bae`
**Scope:** `extern/feasibilityjump/`, `include/kernelpump/`, `src/`, `tests/`, `settings/`, `results/`, the scripts and docs
**Prerequisite:** the `FeasibilityHeuristic` interface (`docs/spec-heuristic-interface.md`, landed)

## 1. Goals and scope

Add Feasibility Jump (FJ) as a third `FeasibilityHeuristic` implementation, registered as `fj`, so it
runs both standalone (`method=fj`) and per bucket inside the Kernel Pump (`kp.heuristic=fj`). Neither
`main.cpp` dispatch nor the KP loop needs to change: the factory provides both roles
(`docs/spec-heuristic-interface.md` §6.5).

There is no `develop` branch. `generic_interface` is the latest development line: it is `main` @
`fc31aa97` plus the pluggable heuristic interface, `method=<name>` dispatch and `submip`.

The research questions this work must answer:

1. **RQ1 Success:** how often does FJ reach a feasible solution, standalone and inside KP, compared with
   FP\*, KP\* and CPLEX?
2. **RQ2 Quality:** how good are those solutions? FJ ignores the objective until it is feasible, so
   quality is the expected weak point.
3. **RQ3 Decomposition:** does KP's kernel and bucket decomposition help FJ converge? Answering this
   needs a control for the restarts that a per-bucket run brings with it (§5.3).

### Why FJ, checked against the `kernel-pump-literature` skill

- The KP paper's future work names **extending the decomposition to other MILP heuristics**. FJ is
  the first LP-free one, so it tests whether the decomposition helps a heuristic that never has a
  fractional point.
- The skill lists FJ as prior art (Luteberget & Sartor, MPC 2023), not as a variant the authors tried.
  Using it is not a novelty claim. **The claim to protect is the decomposition applied to FJ**, which,
  as far as we know, nobody has published.
- KP's mechanism relies on few binaries being active in a feasible solution. FJ starts every variable
  at its lower bound, so the all-zeros scenario is its natural start. That makes KP and FJ a
  plausible fit, but it is a hypothesis for RQ3, not a fact.

### Non-goals

- No change to `FeasibilityPump`, `KernelPump` or `Solution::WriteToFile`. `results/` reads the `.sol`
  format by line offset.
- No change to any committed config. The regression suite must stay at 189/189.
- No multi-threaded FJ. The paper ran FJ on two threads, one on the original model and one on the
  presolved model. The published KP protocol is single-threaded, so `fj` runs on one thread and uses
  one of the two models (§3.4).
- No port of the Rust competition version. Use the C++ reference header only.

## 2. Feasibility Jump in brief

FJ is guided local search over a Lagrangian relaxation. It keeps bounds and integrality, turns every
linear constraint into a weighted violation penalty, and needs no LP solution
([paper](https://doi.org/10.1007/s12532-023-00234-8), `papers/feasibility-jump.pdf`).

### The algorithm (paper, Algorithm 2)

1. Start from any assignment, with every constraint weight `wᵢ = 1`.
2. For each variable, keep one **jump value**: the value in `[lⱼ, uⱼ]` that minimises the weighted
   violation while every other variable is held fixed. Its **score** is how much that move would
   reduce the weighted violation.
3. At each step, sample up to 25 variables with a positive score, pick the best one and make its jump.
   Then update the jump values and scores of its neighbours.
4. At a local minimum (no variable has a positive score), add 1 to the weight of each violated
   constraint and take the best move inside one randomly chosen violated constraint.
5. Stop when every constraint is satisfied, or after a set amount of **effort** without improvement.
   Effort counts inner-loop work, not wall-clock time.

The method is incomplete. It can miss feasible solutions and cannot prove infeasibility, and in the
paper's experiments it ignores the objective entirely. The authors tuned it to run very briefly: 0.6 s
on average.

### What the paper reports

| Item | Value |
| --- | --- |
| Test set | MIPLIB 2017 benchmark, 240 instances |
| Solved by FJ | 123/240, combining both threads |
| Solved on one model only | original model: 6; presolved model: 33. Solved on both: 84 |
| More effort | *"did not see significant improvements"* |
| Competition | 1st place, MIP 2022 Computational Competition |

The presolve split matters here: the presolved model alone accounted for 33 instances. Our
single-threaded `fj` should therefore run on the presolved model by default (§3.4).

### The reference implementation

[SINTEF/feasibilityjump](https://github.com/SINTEF/feasibilityjump) is MIT-licensed; its last commit is
`93f1c2a` (2026-09-25). Only `feasibilityjump.hh` is reusable. `xpress_fj.cc` is an Xpress driver that
relies on global state and detached threads, so it cannot serve as a library.

| Part of `feasibilityjump.hh` | What it means for us |
| --- | --- |
| `FeasibilityJumpSolver(int seed, int verbosity, double weightUpdateDecay)` | one solver object per invocation; seed taken from the global `seed` key |
| `addVar(VarType, lb, ub, objCoeff)` | `Continuous` or `Integer`; binaries are integers on `[0,1]`. Each variable starts at `lb`, with no guard against infinite bounds |
| `addConstraint(RowType, rhs, n, idx, coef, relax_continuous)` | `Equal`, `Lte` and `Gte` only, so range rows must be split. A zero coefficient causes a division by zero |
| `solve(double *initialValues, callback)` | always returns 0 and never stops on its own; the callback must return `Terminate` |
| `FJStatus` | `totalEffort`, `effortSinceLastImprovement`, `numVars`, `solutionObjectiveValue`, `solution` (non-null only when a new feasible point improves the objective) |
| callback timing | on each improving feasible solution, and every 500 000 effort units |
| objective | `objectiveWeight` starts at 0 and rises only at local minima where nothing is violated, so the objective has **no effect before the first feasible point** |
| driver defaults | stagnation limit `maxEffort = 1e8`; stops at the first solution; `weightUpdateDecay = 1.0` (0.9999 with `-e`) |

## 3. Integration design

`fj` is a thin adapter. It reads the model through `MIPModelI` accessors, builds a
`FeasibilityJumpSolver`, and runs it under a wall-clock budget that the callback enforces. It then maps
the result back into the driver's terms. It needs no new backend code, so it behaves the same under
`solver=cpx` and `solver=grb`.

### 3.1 Vendoring

- Copy `feasibilityjump.hh` and `LICENSE` into `extern/feasibilityjump/`, pinned to commit `93f1c2a`.
  Add a `README` that records the commit and every local patch.
- Include the header **only** from `src/fj_heuristic.cpp`. It declares unscoped enums (`Equal`, `Lte`,
  `Integer`, `Continue`, …) and a `FJ_LOG_PREFIX` macro, none of which may leak into other files.
- Before writing the adapter, confirm the header builds as strict C++17 under both MSVC v143 and
  GCC 15 (the `CLAUDE.md` rules).

### 3.2 One small patch to the header

The callback only sees a point when a new **feasible** solution improves the objective. When FJ fails,
KP needs the least-violated assignment instead, for its warm start and its closest-point test. Keep
the patch small and list it in the vendor `README`:

1. Track the assignment with the lowest unweighted total violation, and when it was reached.
2. Expose that assignment and the move count through two read-only accessors.
3. Add an optional total-effort limit alongside the step limit.

If upstream changes make the patch hard to maintain, the fallback is to recompute the violation of
the final assignment only. That loses the best point seen during the run, so it is the second choice.

### 3.3 Building the FJ problem from `MIPModelI`

| Model element | Mapping |
| --- | --- |
| columns | `col()` / `lbs`, `ubs`, `ctypes`, `objcoefs`. `B` and `I` become `Integer`; `C` becomes `Continuous` |
| infinite bounds (±1e20 on CPLEX, `GRB_INFINITY`) | passed as ±∞. The adapter sets the start value (§3.5) rather than relying on `addVar`'s default of `lb` |
| rows `L`, `G`, `E` | `Lte`, `Gte`, `Equal` |
| range rows `R` | split into a `Gte` row and an `Lte` row; a negative range is an error |
| zero coefficients | dropped before `addConstraint`, because the header divides by them |
| objective | built by the same three-way override as FP's initial LP and `submip` (§4.2) |
| `relax_continuous` | `0` by default; exposed as `fj.relaxContinuous` |

FJ needs only the rows, so the adapter is backend-neutral by construction. It must never call
`status()`: there is no MIP solve to classify.

### 3.4 Presolve

`fj` follows `FeasibilityPump::init` exactly. If `mipPresolve=1`, it presolves the model it is given,
works on `presolvedModel()`, and `getSolution()` returns the result of `postsolveSolution()`. The
same code path serves both roles:

- **Standalone:** FJ works on the whole problem, presolved under CPLEX. This matches the paper's
  better thread (+33 instances).
- **Inside KP:** each bucket model has every inactive binary pinned at `ub = 0`, so presolve removes
  those variables. FJ then searches only the active kernel and the current bucket, which is the
  decomposition at work.
- **Gurobi:** `presolve()` does nothing, so `fj` runs on the original model. That is correct
  behaviour, and the experiment report must state it.

### 3.5 Start point and warm start

- **No warm start** (standalone, and KP's first bucket): each variable starts at the value closest to
  0 within its bounds `[l, u]`. For binaries this is the all-zeros scenario. Decided: keep this
  reference (§6.4).
- **KP warm start** (`start.frac` non-empty): use it, but **project it onto the current bounds first**:
  clip each value to the current bounds and round the integer variables. Bounds change between
  buckets, and the header copies the start in without any checks.
- The warm start is FJ's own best point from an earlier bucket, mapped through presolve with
  `presolveSolution()`.

### 3.6 Time, effort and stopping

FJ has no clock, so the callback enforces three stops:

| Stop | Rule |
| --- | --- |
| wall clock | `Terminate` once the elapsed time exceeds `run()`'s `timeLimit`, or on `model->aborted()` (Ctrl-C) |
| stagnation | `Terminate` when `effortSinceLastImprovement > fj.maxEffortNoImpr` (default `1e8`, the driver's value); what happens next is described below the table |
| first solution | `Terminate` on the first feasible point, unless `fj.improve=1` (§4.3) |

The callback fires every 500 000 effort units, which on large instances can be seconds apart. Measure
how far runs overshoot `timeLimit` in testing. If the overshoot is too large, add a shorter callback
interval to the patch in §3.2.

**After a stagnation stop.** A stagnation stop always ends the current FJ run. If
`stopOnNoImprovement` is true (every KP bucket except the last), `run()` returns and KP moves on.
Otherwise (standalone runs, and KP's last bucket), `fj.restarts=1` starts a fresh run with a new seed
and fresh weights, warm-started from the best point so far, and repeats until time runs out. With
`fj.restarts=0` the stagnated run simply continues.

**Time slices.** `fj.timeSlices=k` (standalone only) splits the time limit into `k` equal slices. Each
slice runs a fresh FJ search, warm-started from the best point so far. This reproduces KP's restart
schedule without its decomposition, and it is the control used in §5.3.

### 3.7 Mapping onto the interface

| Member | Implementation |
| --- | --- |
| `readConfig()` | read `fj.*` and the objective keys; called once |
| `init(model)` | presolve as in §3.4 and keep the working model; return `false` only if presolve proves infeasibility |
| `run(tl, start, stopNoImpr)` | build the FJ problem, solve, and map the outcome (table below) |
| `getSolution(x)` | the best feasible point, postsolved, then checked with `isSolutionFeasible(x)` on the caller's model |
| `getPrimalBound()` | that point's value under the **original** objective, not the override |
| `getClosestDist()` | the lowest unweighted total violation seen in this call (`0` when feasible) |
| `getClosestFrac(f)` | the assignment that reached it, postsolved; integral by construction |
| `getIterations()` | moves summed over all calls (an effort total would overflow `int`) |
| `validateConfig()` | rejects `fj.maxEffortNoImpr <= 0` and unsupported key combinations with a visible error, before the silenced handler |

| How the FJ run ended | `foundSolution` | `progressed` | What KP does |
| --- | --- | --- | --- |
| feasible point that passes `isSolutionFeasible` | `true` | `true` | records it and stops |
| stagnation or time limit, no feasible point | `false` | `true` | keeps the point if `dist` beats its running best, then grows the kernel |
| feasible point that fails the row check | `false` | `true` | same as above, and logs a tolerance mismatch |
| empty model, or Ctrl-C | `false` | `false` | absorbs the reference kernel, or breaks out on `aborted()` |

FJ cannot prove infeasibility, so `progressed = false` must never be taken to mean "infeasible".

**State kept across invocations** (the interface requires this to be documented in the header): only
the move total and the restart seed counter carry over. Each call gets a fresh solver, fresh weights
and the driver's time limit. There is no shared budget like `fp.iterLimit`.

### 3.8 Two semantic gaps to state in the results

1. **"Distance" means something different.** For FP, `closest dist` is a fractional point's distance
   to integrality. For FJ, it is an integral point's constraint violation. KP compares distances only
   within one run, so its loop stays correct, but comparing the values across heuristics is
   meaningless.
2. **A run that does not converge reports violation = ∞** (decided 2026-11-23). The violation is used
   only inside the run, for KP's warm start and its closest-point test, and is never reported as
   partial credit. There is no new `.sol` line and no change to `results/`. On FJ's integral closest
   point, the `real integrality gap` and `projection integrality gap` fields would read 0, so the
   results treat them as not applicable for `fj` runs. A failed `fj` run counts as a failure.

## 4. Configuration

Registering `fj` provides both entry points without new dispatch code: `method=fj` runs it standalone
and `kp.heuristic=fj` runs it per bucket. Every committed config keeps `kp.heuristic=feaspump`.

### 4.1 New `fj.*` keys

| Key | Default | Meaning |
| --- | --- | --- |
| `fj.maxEffortNoImpr` | `1e8` | effort without improvement before a stagnation stop (the driver's value); what happens next depends on `stopOnNoImprovement` (§3.6) |
| `fj.maxEffort` | `0` | cap on total effort per call; `0` means no cap, so the wall clock decides |
| `fj.weightUpdateDecay` | `1.0` | `1.0` gives additive weights, as in the paper; `0.9999` is the driver's `-e` variant |
| `fj.restarts` | `1` | after a stagnation stop with time left, start again with a new seed and fresh weights |
| `fj.timeSlices` | `1` | standalone only: split the time limit into this many fresh, warm-started runs (§3.6) |
| `fj.improve` | `0` | `0` stops at the first feasible point; `1` keeps improving the objective until stagnation or the time limit (§4.3) |
| `fj.relaxContinuous` | `0` | the header's `relax_continuous` flag |
| `fj.verbosity` | `0` | the header's log level; forced to 0 under `-DSILENT_EXEC` |

`fj` also honours the global keys `seed`, `timeLimit` and `mipPresolve`. An unknown `fj.*` key is
silently ignored, as everywhere else in this codebase.

### 4.2 The objective: no new key

Follow `submip` §6.1, so that KP⁻, KP⁺ and KP⁺ᐟ⁻ mean the same whichever heuristic runs inside them.
FJ's objective coefficients come from the same chain as FP's initial LP, applied in the same order:

| Key | FJ objective |
| --- | --- |
| `fp.forceNullObjectiveInitialLP=1` | all zeros |
| else `fp.forceSumVarsObjectiveInitialLP=1` | 1 on every binary; `...MaxSense=1` switches to maximisation |
| else | the original objective |
| `fp.reverseObjectiveFunction=1` | applied last: flip the sense of the result |

Write it as the same `if / else if` chain that FP uses. In the `submip` work, merging these branches
introduced a precedence bug.

### 4.3 What the objective can and cannot do in FJ

In FJ, the objective has **no effect before the first feasible point** (§2). With `fj.improve=0`, the
override chain therefore changes nothing in the search. KP⁻ and KP⁺ still differ under `fj`, but only
through the kernel-building LP (`kp.*` keys).

That makes `fj.improve` the main lever for RQ2. In a standalone run, `fj.improve=1` spends the time
left after the first solution on improving it. Inside KP, the loop stops at the first solution the
heuristic returns, so `fj.improve=1` acts only inside the bucket that found it. `main.cpp` then
reoptimizes every solution by fixing the integer variables and re-solving the LP, which repairs the
continuous part for every method alike. Decided: both modes are tested and reported (§6.4).

### 4.4 New config files

Each new file copies its base and changes only the lines listed. None of them goes into
`settings/settings_list.txt`, so the paper's reproduction script is unaffected.

| Config | Base | Changes | Role |
| --- | --- | --- | --- |
| `config29` | `config1` (FP\*) | `method=fj` | **FJ**: standalone, first solution |
| `config30` | `config29` | `fj.improve=1` | **FJ-imp**: standalone, improving |
| `config31` | `config16` (KP\*) | `kp.heuristic=fj` | **KP\*+FJ** |
| `config32` | `config17` (KP⁻) | `kp.heuristic=fj` | **KP⁻+FJ** |
| `config33` | `config21` (KP⁺) | `kp.heuristic=fj` | **KP⁺+FJ** |
| `config34` | `config29` | `fj.timeSlices=N`, with N taken per instance and seed from the KP\*+FJ run (§5.3) | **FJ-restart**: the control for RQ3 |
| `config35` | `config31` | `fj.improve=1` | **KP\*+FJ-imp** |
| `config36` | `config32` | `fj.improve=1` | **KP⁻+FJ-imp** |
| `config37` | `config33` | `fj.improve=1` | **KP⁺+FJ-imp** |

## 5. Experimental design

All methods run on the same machine with the same CPLEX version, single-threaded, so that FJ times
are comparable with the baselines. The published tables come from different hardware, so FP\*, KP\*,
KP⁻, KP⁺ and CPLEX_feas are **re-run**, not copied from the paper.

### 5.1 Methods and the question each comparison answers

| Comparison | Configs | Answers |
| --- | --- | --- |
| FJ vs FP\* vs CPLEX_feas | 29, 1, 18 | RQ1: is FJ alone competitive on our test sets? |
| KP\*+FJ vs KP\* | 31, 16 | RQ1: does FJ beat FP as KP's inner heuristic? |
| KP⁻+FJ, KP⁺+FJ vs KP⁻, KP⁺ | 32, 33, 17, 21 | RQ1: does that hold for the objective extremes? |
| each FJ method in both modes: first solution vs improving | 29/30, 31/35, 32/36, 33/37 | RQ2: how much does improvement mode gain, and what does it cost in success and time? Every table reports both modes |
| KP\*+FJ vs FJ-restart vs FJ | 31, 34, 29 | RQ3: decomposition vs restarts vs neither |

### 5.2 Metrics

| Metric | Source | Used for |
| --- | --- | --- |
| success: a feasible solution within the limit; otherwise violation = ∞ | `.sol` `status` | RQ1, RQ3 |
| time to first feasible solution | `.sol` `total time` | RQ1, RQ3 |
| primal gap, `abs(z − z*) / max(abs(z), abs(z*))` (0 if both are 0), on the reoptimized value | `.sol` + `results/bestKnownBounds.csv` | RQ2 |
| FJ moves to the first feasible solution | `# iterations` | RQ3 |
| last bucket visited, and active binaries / all binaries at the end | `.sol` | RQ3 |
| binaries at 1 in the solution / all binaries | `.sol` | RQ3: groups instances by the structure KP exploits |

The gap is reported only over the instances that **every compared method solved**, together with how
many instances that is. Otherwise a method that solves only easy instances would look better.

### 5.3 RQ3: separating decomposition from restarts

KP+FJ differs from FJ in two ways at once. It **decomposes** the problem into small sub-problems,
ordered by the LP, with the kernel growing as it goes. It also **restarts** FJ at every bucket, with
fresh weights and a time slice of `δ = T/(|Λ|+1)`. A gain could come from either one.

The FJ-restart control (`config34`) keeps the restarts and drops the decomposition. It runs
standalone FJ with `fj.timeSlices=N`, where `N` is the number of FJ calls KP\*+FJ made on that
instance and seed. `script_run_experiments.sh` reads `N` from the KP run's `.sol` file and passes it
as an inline override.

| Result | Reading |
| --- | --- |
| KP\*+FJ > FJ-restart ≈ FJ | the decomposition helps FJ converge |
| KP\*+FJ ≈ FJ-restart > FJ | the gain comes from restarts; KP adds nothing |
| KP\*+FJ > FJ-restart > FJ | both contribute |
| KP\*+FJ ≤ FJ | the decomposition hurts FJ, for example because the slices are too short |

If the decomposition helps, the effect should be largest where few binaries are active in the
solution. Report RQ3 per test set and by that share (§5.2): R-STOP-DP and TTP-PV are where KP's wins
were statistically significant.

A second control would keep KP's bucket sizes but fill the buckets at random (a new
`kp.shuffleBuckets` key). It would separate the LP-guided order from the effect of growing the kernel
at all. Because it changes KP code, it is a follow-up and not part of this work (decided, §6.4).

### 5.4 Stages

| Stage | Instances | Limit | Seeds | Purpose |
| --- | --- | --- | --- | --- |
| E0 smoke | ~10 from the golden set | 60 s | 1 | check the adapter in both roles, with no row-check failures |
| E1 pilot | all three sets | 300 s | 1 | first answers to RQ1–RQ3; settle `fj.maxEffortNoImpr` and `fj.improve` |
| E2 full | MIPLIB 2017 (233), R-STOP-DP (48), TTP-PV (36) | 3600 s | 3 | the published protocol; the numbers for a paper |

The FJ paper reports little gain from more effort, so the 300 s pilot may already show the pattern.
E2 is still needed to compare with the published KP results.

### 5.5 Statistics

- **Success:** Cochran's Q across all methods, then pairwise McNemar tests with Holm correction at 5%,
  as in the paper (`results/src/McNemar_test.py`).
- **Time and gap:** pairwise Wilcoxon signed-rank tests on the instances both methods solved, with
  Holm correction. This test does not exist yet; add it next to the McNemar script.
- **Time profiles:** performance profiles of time to first feasible solution, as in the FJ paper.
- Count each seed as a separate observation for McNemar, as the paper does, and say so in the
  write-up.

## 6. Plan, risks and decisions

### 6.1 Milestones

| # | Content | Done when |
| --- | --- | --- |
| **F0** | vendor `feasibilityjump.hh` at `93f1c2a` + `LICENSE`; apply the §3.2 patch | it builds as strict C++17 on MSVC v143 and GCC 15 |
| **F1** | `fj_heuristic.{h,cpp}`, factory registration, problem build (§3.3), presolve (§3.4), stops (§3.6) | `method=fj` solves a toy model and a golden-set instance |
| **F2** | KP role: warm-start projection, closest point and distance, restarts and time slices | `kp.heuristic=fj` runs end to end on ~5 instances |
| **F3** | contract tests (§6.2) | all pass, alongside the existing 79 checks |
| **F4** | regression suite | 189/189 unchanged; nothing on the default path changes |
| **F5** | configs 29–37, `fj.timeSlices` support in `script_run_experiments.sh`, Wilcoxon script, docs (`TESTSET.md`, `CLAUDE.md`, the three skills) | E0 runs from the scripts |
| **F6** | E0 smoke, then E1 pilot (§5.4) | first answers to RQ1–RQ3; `fj.*` defaults settled |
| **F7** | E2 full protocol | results tables and statistical tests for the paper |
| **F8** | Gurobi check on toy models | needs the academic licence; never the work licence |

F0–F7 run on a machine with a CPLEX licence.

### 6.2 Tests to add to `tests/heuristic_interface_test.cpp`

1. `fj` resolves case-insensitively and works for both `method=` and `kp.heuristic=`.
2. A feasible toy MIP gives `{true, true}`, and the returned point passes `isSolutionFeasible`.
3. An infeasible toy MIP gives `{false, true}` at the time limit, never `{false, false}`: FJ cannot
   prove infeasibility.
4. Problem build: a range row becomes two rows, a zero coefficient is dropped, and a variable with an
   infinite bound gets a finite start.
5. Warm start: a start outside the current bounds is clipped and rounded before `solve()`.
6. Time: on a model where callbacks are far apart, the overshoot past `timeLimit` stays below a
   chosen bound.
7. Determinism: the same seed gives the same moves and the same point; a different seed gives a
   different run.
8. Objective chain: the precedence is checked directly in the code path, as for `submip`, because
   that is where the earlier bug was.
9. Presolve round trip: the closest point survives `postsolveSolution` and `presolveSolution` with the
   right size.
10. Stagnation: with `stopOnNoImprovement=true` the run returns; with `false` and `fj.restarts=1` it
    restarts with a new seed.

### 6.3 Risks

| Risk | Mitigation |
| --- | --- |
| with callbacks only every 500 000 effort units, runs on large models overshoot the time limit | test 6; shorten the interval in the patch if needed |
| FJ's 1e-5 tolerances disagree with `isSolutionFeasible` | report success only after the row check, and log each mismatch |
| big-M rows and infinite bounds make jump values extreme | use a finite start (§3.5); count runs that start on an infinite bound |
| short KP slices stop FJ before it ever stagnates | that is a real result for RQ3; `last bucket visited` shows it |
| the header patch drifts from upstream | keep it small, list it in the vendor `README`, pin the commit |
| results across heuristics misread through `closest dist` or the gap fields | §3.8, stated in every results table |

### 6.4 Decisions (2026-11-23)

| # | Question | Decision |
| --- | --- | --- |
| 1 | Start point of FJ in KP's first bucket | keep the all-zeros reference (§3.5); no KP change |
| 2 | `fj.improve` default | test and report both modes for every FJ method (configs 29–37) |
| 3 | Violation in the `.sol` file | none; a run that does not converge reports violation = ∞ (§3.8) |
| 4 | Shuffled-bucket control | follow-up, outside this work (§5.3) |
| 5 | Config numbers | 29–37, as listed in §4.4 |

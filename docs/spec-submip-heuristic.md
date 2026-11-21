# Spec: a sub-MIP feasibility heuristic (`submip`), portable across CPLEX and Gurobi

**Status:** proposed
**Branch:** to be created from `generic_interface`
**Scope:** `include/kernelpump/`, `src/`, `tests/`, `settings/`, docs
**Prerequisite:** the `FeasibilityHeuristic` interface (`docs/spec-heuristic-interface.md`, M0–M6, landed)

## 1. Goal

Add a second `FeasibilityHeuristic` implementation, registered as `submip`, that looks for integer
feasibility by handing the restricted problem to the backend's own MIP solver and stopping at the
first feasible solution. It must be usable in both roles the interface now supports:

- `kp.heuristic=submip` — driven per bucket by the Kernel Pump, so kernel growth is guided by an
  exact sub-solver instead of a pump;
- `method=submip` — standalone over the whole problem.

**Hard requirement from the outset: it must run under both `solver=cpx` and `solver=grb`.** The two
backends differ in ways that silently break naive code (Section 4); designing for that now is much
cheaper than retrofitting it.

### Why this heuristic, from the literature

Checked against the `kernel-pump-literature` skill:

- The KP paper's own future work lists **"extending decomposition to other MILP heuristics"**. The
  interface exists for that; this is the first heuristic that exercises it.
- It is also the natural second member of a **portfolio**, which the paper lists as a future
  direction and which Salvagnin, Roberti & Fischetti (MPC 2024) found far more robust than any
  single strategy on a heterogeneous testbed.
- It is **not** one of the variants the authors tried and rejected (those were all about the
  objective driving the kernel, or the ordering criterion — see the skill's rejected-ideas table).
- A sub-MIP node limit of **500** is the precedent set by Mexi, Berthold & Salvagnin (EJCO 2023) for
  sub-MIP stages inside FP, matching SCIP conventions. Use it as the default rather than inventing
  one.

It also gives the paper a baseline it currently lacks: KP's decomposition measured against an exact
sub-solver, separating "the decomposition helps" from "the pump helps".

## 2. Non-goals

- No change to `FeasibilityPump`, to `KernelPump`'s loop, or to the `FeasibilityHeuristic`
  interface itself, except the one optional extension in §7 (MIP starts), which is a later milestone.
- No change to the `.sol` format (`results/` parses it by line offset — see
  `docs/spec-heuristic-interface.md` §7.4).
- No new default behaviour: every committed config keeps `kp.heuristic=feaspump`, so the regression
  suite must stay at 190/190.
- Not a claim of research merit. This spec delivers a correct, portable implementation and the
  experiments that would support or refute it; whether it earns a place in a paper is a separate
  question.

## 3. The algorithm

Per invocation, on the model the driver supplies (for KP: the current kernel ∪ bucket, with every
inactive binary pinned by `ub = 0`):

1. **`init(model)`** — clone the model and keep the clone. Return `true`.
   No presolve (see §4.2), so no crush/uncrush and no vector mapping anywhere.
2. **`run(timeLimit, start, stopOnNoImprovement)`** — on the clone:
   - `Threads = 1` unless `multiThreading`;
   - `seed` from the global `seed` key;
   - `DblParam::TimeLimit = timeLimit`;
   - `IntParam::SolutionLimit = 1` — stop at the first feasible solution;
   - `IntParam::NodeLimit = submip.nodeLimit` (default 500; 0 means unlimited);
   - `IntParam::Emphasis = 1` if `submip.feasEmphasis` (default on);
   - **replace the objective exactly as FP's initial LP does** (§6.1);
   - `mipopt()`;
   - classify the outcome **by the portable predicates only** (§4.1), never by `status()`.
3. **Outcome mapping** onto `HeuristicOutcome`:

| Sub-MIP ended | `foundSolution` | `progressed` | What KP does with it |
| --- | --- | --- | --- |
| a feasible point is available | `true` | `true` | records the solution, absorbs the reference kernel, stops |
| proven infeasible | `false` | `false` | absorbs the whole reference kernel permanently — the correct reaction |
| limit hit (time/nodes), no solution | `false` | `true` | continues to the next bucket |
| Ctrl-C (`aborted()`) | `false` | `false` | KP's loop breaks on `model_->aborted()` anyway |

`stopOnNoImprovement` is ignored: there is no stagnation notion in a branch-and-bound call. Document
that, because KP passes `false` only for the last bucket and an implementer may expect it to matter.

4. **Double-check before reporting success.** `isPrimalFeas()` alone is not enough (§4.1); verify the
   retrieved point with `model->isSolutionFeasible(x)`, which is implemented once in `mipmodel.h`
   against the row data and is therefore backend-independent. Report success only if both hold. This
   mirrors `main.cpp`'s own final feasibility re-check, which is the project's correctness habit.

## 4. The portability contract — CPLEX and Gurobi

This section is the reason the spec exists. Everything below was verified by reading
`src/cpxmodel.cpp` and `src/grbmodel.cpp`.

### 4.1 `mipopt()`'s return value is meaningless — on both backends

```cpp
bool CPXModel::mipopt() { ... auto solve_status = this->status(); return true; }   // src/cpxmodel.cpp
bool GRBModel::mipopt() { GRBoptimize(model); dirty = false; return true; }        // src/grbmodel.cpp
```

Both ignore the solve outcome and return `true` unconditionally. **Never branch on it.**

`status()` is equally unusable: it returns the raw solver code, `CPX_STAT_*` versus `GRB_*`, which do
not share values or meanings. Branching on `status()` is the single easiest way to write code that
"works" on CPLEX and misbehaves on Gurobi.

The portable predicates, both implemented with deliberately aligned semantics:

| Call | Meaning | Note |
| --- | --- | --- |
| `isInfeasibleOrTimeReached()` | infeasible, or stopped by a limit with nothing usable | `GRBModel`'s implementation explicitly mirrors CPLEX's rule: for a MIP, a time-out counts only when `SolCount == 0` |
| `isPrimalFeas()` | a usable point exists | Gurobi checks `SolCount > 0` and, if stopped early, that `ConstrVio`/`BoundVio` are within `FeasibilityTolerance` |
| `aborted()` | Ctrl-C | both install a handler |
| `isSolutionFeasible(x)` | row-by-row check | non-virtual helper in `mipmodel.h`, identical on every backend |

`main.cpp`'s existing `method=solver` path is the template: it uses
`foundSolution = !(model->isInfeasibleOrTimeReached())` and nothing else. Follow it, but tighten it
with `isPrimalFeas()` **and** `isSolutionFeasible()` as in §3.4 — the existing one-predicate test
would call "found" a solve that ended for some other reason with no solution, and a heuristic that
reports a phantom solution to KP would corrupt the run.

**Distinguishing "proven infeasible" from "limit reached, nothing found"** needs care, because
`isInfeasibleOrTimeReached()` deliberately conflates them. Portable approach: treat
`isInfeasibleOrTimeReached() && !isPrimalFeas()` as "no solution", and separate the two cases by
whether a limit was actually in force and reached — i.e. compare elapsed time against `timeLimit`
and the node count (`IntAttr::Nodes`) against `submip.nodeLimit`. If neither limit was reached and
no solution exists, the sub-MIP was solved to completion and is infeasible → `progressed = false`.

That inference is slightly conservative: a sub-MIP that is genuinely infeasible *and* hit a limit is
reported as "limit reached", so KP continues rather than absorbing the bucket. That is the safe
direction — it never absorbs on a guess. **If this proves too coarse, the clean fix is a new
`MIPModelI` predicate (`isProvenInfeasible()`), implemented as `CPX_STAT_INFEASIBLE`/
`CPXMIP_INFEASIBLE` and `GRB_INFEASIBLE`, which keeps the solver codes behind the abstraction where
they belong.** Prefer that over leaking status codes into the heuristic.

### 4.2 Gurobi has no MIP presolve — so do not presolve

`GRBModel::presolve()` warns once and returns `true` without doing anything; `presolvedModel()`
yields nothing, and `postsolveSolution`/`presolveSolution` are identity. The C API has no
crush/uncrush.

Consequences, all of which make the design *simpler*:

- `submip` must not presolve. It works on the model as handed over, so solution vectors are already
  in the driver's space and need no mapping — exactly the class of bug §7.3 of the interface spec
  warns about.
- Under `solver=grb`, KP as a whole already runs unpresolved: `KernelPump::Init` sees
  `presolvedModel() == nullptr`, logs "no reductions", clones, and leaves `has_presolve_` false.
  That path is already correct; `submip` inherits it.
- `FeasibilityPump` relaxes only its **own** copy (`model = premodel;` then `switchToLP()`), never
  the model the driver passed, so after an FP bucket the driver's model still carries its integer
  types. A sub-MIP heuristic can therefore follow an FP one in the same run.

### 4.3 Every parameter and attribute needed is supported on both

Verified against `GRBModel::intParam`, `dblParam`, `intAttr` and `seed`:

| Need | `MIPModelI` | CPLEX | Gurobi |
| --- | --- | --- | --- |
| stop at first solution | `IntParam::SolutionLimit` | `CPX_PARAM_INTSOLLIM` | `SolutionLimit` |
| node cap | `IntParam::NodeLimit` | node limit | `NodeLimit` (as a double; `INT_MAX` → `GRB_INFINITY`) |
| time cap | `DblParam::TimeLimit` | ✓ | ✓ |
| feasibility emphasis | `IntParam::Emphasis` | emphasis = 1 | `MIPFocus = 1` |
| single thread | `IntParam::Threads` | ✓ | ✓ |
| reproducibility | `seed(int)` | random seed | `Seed` |
| work done | `IntAttr::Nodes` | node count | `NodeCount` |
| relax to LP (only if §7.2 is built) | `switchToLP()` | ✓ | ✓ |

So the heuristic needs **no new backend code** for the core milestones. That is the main design
constraint discharged.

### 4.4 Gurobi cannot be compiled or tested on the current dev machine

`find_package(GUROBI)` fails here: `GUROBI_INCLUDE_DIR-NOTFOUND`, and there is no `gurobi_c.h`
anywhere on disk. The conda/pip `gurobipy` package ships `libgurobi130.so` but no headers, so the
full distribution tarball plus `-DGUROBI_HOME=...` is required.

Also note the only non-work licence available locally is the bundled pip one, capped at
**2000 × 2000**, which cannot run MIPLIB instances. So:

- the Gurobi path must be written to the contract above and **reviewed**, not assumed tested;
- it can be exercised on toy models within the size-limited licence, which is exactly what
  `tests/grbmodel_test.cpp` already does — follow that pattern (explicit licence file, tiny models);
- any claim that `solver=grb` works end-to-end on real instances must come from a machine with a
  full licence. State that in the final report rather than implying coverage.

## 5. Mapping onto the interface

| Member | Implementation |
| --- | --- |
| `readConfig()` | read the `submip.*` keys (§6) |
| `init(model, ctype)` | clone, store, `true` |
| `run(tl, start, _)` | §3 |
| `foundSolution()` | the last outcome's `foundSolution` |
| `getSolution(x)` | the stored point, already in the driver's space |
| `getPrimalBound()` | `objval()` when a solution exists, else `±INFBOUND` by `objSense()` |
| `getClosestDist()` | pass through `start.dist` (§7.1) |
| `getClosestFrac(f)` | pass through `start.frac` (§7.1) |
| `getIterations()` | accumulated `IntAttr::Nodes` across invocations |
| `name()` | `"submip"` |

**Cross-invocation state** (the interface requires this to be documented in the header): the node
total behind `getIterations()` accumulates; nothing else carries over, and in particular there is no
shared budget analogous to FP's `fp.iterLimit`. Each invocation gets the driver's time limit and a
fresh node limit. Say so explicitly, because a reader who knows FP will assume otherwise.

## 6. Configuration

| Key | Default | Meaning |
| --- | --- | --- |
| `submip.nodeLimit` | `500` | node cap per invocation; `0` = unlimited. Default follows Mexi et al. (EJCO 2023) |
| `submip.feasEmphasis` | `1` | set `IntParam::Emphasis` (CPLEX feasibility emphasis / Gurobi `MIPFocus`) |
| `submip.reportRelaxationPoint` | `0` | §7.2, later milestone |
| `submip.useMipStart` | `0` | §7.3, later milestone |

No existing config is touched. An unknown `submip.*` key is silently ignored
(`FileConfig::getAsString` returns the default for a missing entry), as everywhere else in this
codebase — worth knowing when a setting appears to have no effect.

### 6.1 The sub-MIP objective is **not** a new key

`submip` must **not** introduce an objective switch of its own. It reuses the mechanism the KP
variants already use, so that KP⁻, KP⁺ and KP⁺ᐟ⁻ mean the same thing whichever inner heuristic is
driving — the variant is a property of the configuration, not of the heuristic.

The sub-MIP solve is the `submip` counterpart of FP's **initial LP**, so it applies the same
three-way override that `FeasibilityPump::solveInitialLP` applies (`src/feaspump.cpp`), in the same
precedence order, to its own clone:

| Key | Effect on the sub-MIP objective |
| --- | --- |
| `fp.forceNullObjectiveInitialLP=1` | all objective coefficients 0, offset 0 — pure feasibility |
| `fp.forceSumVarsObjectiveInitialLP=1` | coefficient 1 on every binary, 0 elsewhere, offset 0; sense `MIN` unless `fp.forceSumVarsObjectiveInitialLPMaxSense=1`, then `MAX` |
| neither | the original MILP objective, untouched |
| `fp.reverseObjectiveFunction=1` | afterwards, flip whatever sense is in force |

Precedence is null → sum-of-binaries → original, exactly as in FP; `reverseObjectiveFunction` is
applied last, on top of the result.

Why this is the right objective family rather than an arbitrary one: minimising `Σ xᵢ` over the
binaries **is** minimising the L1 distance to the all-zeros scenario, and maximising it is
minimising the distance to all-ones, since `Σ(1 − xᵢ) = |B| − Σ xᵢ`. So these keys select which
extreme binary scenario the sub-problem is pulled towards — the same intuition the paper builds KP⁻
and KP⁺ on, and the reason the kernel-building LP has its own mirror of the same switch
(`kp.buildKernelBasedOnSumVarsObjective`, `kp.buildKernelBasedOnSumVarsObjectiveMaxSense`,
`kp.buildKernelBasedOnNullObjective`, `kp.reverseObjectiveFunction`).

Consequence worth stating: with `config17`/`config21`/`config23` and nothing but
`kp.heuristic=submip` added, the run keeps the published variant's semantics on both the
kernel-building LP (`kp.*`, unchanged) and the sub-problem (`fp.*`, now honoured by `submip`).

Two notes for the implementer:

- These keys are spelled `fp.*` but now govern a non-FP heuristic. Read them under those names in
  v1 — renaming would silently change what existing configs mean. If the misnomer becomes annoying,
  add neutral aliases (`heur.*`) that **fall back** to the `fp.*` names, never a replacement.
- FP does not restore the original objective after its initial LP, and neither should `submip` —
  there is only one solve per invocation, so the question does not arise, but do not "fix" it in FP
  while looking.

The paper's warning still applies and is now a reporting obligation rather than a key: the authors
found that discarding the original objective holds the success rate but *"substantially degrades
solution quality"*. Any experiment run with `fp.forceNullObjectiveInitialLP=1` must report solution
quality, not just success.

## 7. Three deliberate limitations, and how to lift them

### 7.1 It reports no closest point of its own

A branch-and-bound call that fails produces no "closest to feasible" fractional point. So
`getClosestFrac`/`getClosestDist` **pass the driver's own warm start straight back**. This is honest
(the heuristic did not improve on it) and harmless: KP's improvement test is
`lessThan(curr_dist, closest_dist_)`, which is false for an unchanged value, so the kernel grows
only via the bucket-absorption branch.

Two consequences to state plainly:

- With `kp.forceBucketVarsIntoKernel=0` (not the published setting) KP would grow the kernel *only*
  on closest-point improvements — which `submip` never reports — so the kernel would never grow and
  every bucket would be tried against the initial kernel alone. **`submip` is only sensible with
  `kp.forceBucketVarsIntoKernel=1`.** The implementation should refuse the combination loudly rather
  than run a silently pointless configuration.
- On a KP run that finds nothing, `closest_frac_` stays empty, so the `.sol` `real integrality gap`
  is computed on an empty vector (`computeIntegralityGap` then charges 1.0 per integer variable).
  Defined, but not comparable with FP's gaps. Note it where results are interpreted.

### 7.2 Optional: report the relaxation point (`submip.reportRelaxationPoint`)

To make the gap fields meaningful, solve the LP relaxation of the sub-MIP on a second clone via
`switchToLP()` (implemented on both backends) and report that point with its integrality distance
from `computeIntegralityGap`. Costs one extra LP per invocation. Later milestone, off by default.

### 7.3 Optional: seed the sub-MIP with the warm start (`submip.useMipStart`)

`MIPModelI` has **no MIP-start facility**, so this needs one new interface method, e.g.
`virtual void mipStart(const std::vector<double> &x) = 0;`, implemented as:

- CPLEX: `CPXaddmipstarts` with effort level `CPX_MIPSTART_AUTO`;
- Gurobi: `GRBsetdblattrarray(model, "Start", 0, n, values)`.

The warm start is fractional, so supply it as a **partial** start: pass only components already
integral and mark the rest unknown (`GRB_UNDEFINED` for Gurobi; omit the index for CPLEX). Adding a
pure-virtual to `MIPModelI` forces every backend — including the non-production XPRESS/SCIP/PDLP
classes — to implement or stub it, so prefer a virtual with a no-op default. Later milestone.

## 8. Files

| File | Change |
| --- | --- |
| `include/kernelpump/submip_heuristic.h` | **new** — the class and its documented cross-invocation state |
| `src/submip_heuristic.cpp` | **new** — implementation + factory registration (`SUBMIP_FACTORY_RECORDER`, the idiom from `src/feaspump_heuristic.cpp`) |
| `CMakeLists.txt`, `results/CMakeLists.txt` | add the new TU to `libkp` |
| `tests/heuristic_interface_test.cpp` | extend (§9) |
| `settings/all/config26.cfg` *(new number)* | a KP+submip config for experiments; existing configs untouched |
| `TESTSET.md`, `CLAUDE.md`, the three skills | document the new heuristic and keys once it lands |
| `src/feaspump.cpp`, `include/kernelpump/feaspump.h`, `src/kernelpump.cpp` | **none** |

## 9. Milestones

| # | Content | Done when |
| --- | --- | --- |
| **S1** | `submip` class + registration, CPLEX path only | `method=submip` and `kp.heuristic=submip` both run; regression suite still 190/190 |
| **S2** | portability review against §4, plus the `isProvenInfeasible()` predicate if §4.1's inference proves too coarse | outcome classification uses no solver status code anywhere |
| **S3** | contract tests (§10) | all pass |
| **S4** | Gurobi verification on a machine with a full licence | `solver=grb` with `kp.heuristic=submip` runs on real instances |
| **S5** *(optional)* | `submip.reportRelaxationPoint` | gap fields meaningful for submip runs |
| **S6** *(optional)* | `mipStart` interface method + both backends | measurable effect, or dropped |
| **S7** | **baseline experiment** (in scope): KP+submip vs KP+feaspump on the same configs, plus FP\* and the plain solver for reference | a defensible answer to "decomposition or pump?", reporting solution quality as well as success rate |
| **S8** *(follow-up, not part of this work)* | **contender experiments**: submip as a method in its own right and as a portfolio member, on the full published protocol | deferred — see §12.3 |

The implementation work ends at S7. S8 is a separate piece of research, listed only so the scope
boundary is explicit: do not hold S1–S7 open waiting for it.

## 10. Test plan

Extend `tests/heuristic_interface_test.cpp` (CPLEX-gated, currently 53 checks):

1. **Factory and dispatch.** `submip` resolves via `resolveHeuristicName`, case-insensitively, and
   appears in `registeredHeuristicNames()`. Both `method=submip` and `kp.heuristic=submip` are
   accepted; this is free, given the M6 machinery.
2. **Outcome mapping**, on tiny purpose-built models where the answer is known by hand:
   - a feasible toy MIP → `{true, true}`, and the returned point passes `isSolutionFeasible`;
   - an infeasible toy MIP (e.g. `x1 + x2 = 5` over two binaries) → `{false, false}`;
   - `submip.nodeLimit=1` on something that needs search → `{false, true}`, not `{false, false}`.
   That third case is the one that catches a status-code regression, and it is exactly the case that
   would behave differently on Gurobi if `status()` were used.
3. **Lifecycle**, reusing the existing `RecordingHeuristic` harness around `submip`: one instance per
   KP run, `readConfig()` once, `init`+`run` per bucket.
4. **No shared budget.** Node counts per invocation are each bounded by `submip.nodeLimit`, in
   contrast to FP's globally shared `fp.iterLimit` — assert the contrast, since it is the single
   most likely wrong assumption.
5. **The refusal** of `kp.heuristic=submip` with `kp.forceBucketVarsIntoKernel=0` (§7.1).
6. **Standalone equivalence.** With `submip.nodeLimit=0` and `submip.feasEmphasis` matched to
   `mipFeasEmphasis`, `method=submip` should find a solution on the same toy instances as
   `method=solver`. Not a bitwise check — both stop at the first feasible solution found by the same
   solver, so statuses should agree even when the points differ.
7. **The objective override is applied.** On a toy model whose original objective and whose
   `Σ binaries` objective prefer *different* first feasible solutions, check that
   `fp.forceSumVarsObjectiveInitialLP=1` with `MaxSense=0` yields the all-zeros-leaning solution and
   `MaxSense=1` the all-ones-leaning one, and that with neither key the original objective's
   solution comes back. This is what keeps KP⁻/KP⁺ meaningful under `submip`.
8. **Gurobi toy tests** (S4, Gurobi-gated): the outcome-mapping cases from (2) re-run under
   `solver=grb` with an explicit licence file, following `tests/grbmodel_test.cpp`. All fit the
   2000 × 2000 size-limited licence.
9. **Regression suite unchanged**: `kernel-pump-regression` stays 190/190, since nothing on the
   default path moves.

## 11. Risks

| Risk | Mitigation |
| --- | --- |
| outcome classified via `status()`, working on CPLEX and misbehaving on Gurobi | §4.1; test case (2c); S2 review gate |
| "found a solution" reported on a solve that found none | require `isPrimalFeas()` **and** `isSolutionFeasible()` |
| a presolve call sneaks in and breaks under Gurobi, or needs crush/uncrush | §4.2: never presolve; work on the model as given |
| parameters leak into the driver's model across buckets | operate on a clone; never touch the caller's object |
| silently useless with `forceBucketVarsIntoKernel=0` | refuse loudly (§7.1) |
| Gurobi path shipped untested and assumed working | §4.4: say what was and was not exercised |
| sub-MIP per bucket is far slower than a pump, so KP never reaches later buckets | that is a legitimate negative result; `submip.nodeLimit` is the control, and `last bucket visited` in the `.sol` measures it directly |

## 12. Decisions from the author

Answered 2026-11-21; recorded here so the implementer does not re-open them.

1. **Objective inside the sub-MIP — settled: no new key.** It follows the objective switch the KP
   variants already use, i.e. the same `fp.forceNullObjectiveInitialLP` /
   `fp.forceSumVarsObjectiveInitialLP` / `fp.forceSumVarsObjectiveInitialLPMaxSense` /
   `fp.reverseObjectiveFunction` triple that `FeasibilityPump::solveInitialLP` applies, with the same
   precedence. Minimising `Σ xᵢ` over the binaries is the L1 distance to the all-zeros scenario and
   maximising it the distance to all-ones, so these keys already express "pull the sub-problem
   towards one of the two extreme binary scenarios". Full detail in §6.1. The invented
   `submip.nullObjective` key is **dropped**.

2. **Node limit — settled: 500**, following Mexi, Berthold & Salvagnin (EJCO 2023). Exposed as
   `submip.nodeLimit` so it can be swept later, but 500 is the default and needs no justification
   beyond the citation.

3. **Experimental role — settled: both, but only the baseline is in scope here.** Nothing in the
   implementation depends on this; it decides which experiments run when. The two readings:

   - **As a baseline.** Run KP with `submip` and compare against KP with `feaspump`. Both use the
     same kernel and the same buckets, so the only thing that differs is what solves each
     sub-problem. That isolates how much of KP's reported gain comes from *the decomposition* versus
     from *the pump*: if KP+submip already beats plain FP, the decomposition is carrying the result;
     if it does not, the pump is essential to it. This answers a reviewer-style question about the
     published method.
   - **As a contender.** Treat KP+submip as another method in the comparison table alongside FP\*,
     KP\*, KP⁻, KP⁺ and CPLEX, and as a member of the parallel portfolio the paper's future work
     proposes — the fix-propagate-repair paper's finding was that a portfolio is far more robust than
     any single strategy, and KP⁻/KP⁺/KP\* are already a natural portfolio that `submip` would
     extend.

   **Both are wanted, and they are sequenced.** The baseline comparison is part of delivering
   `submip` (S7): it is cheap — KP+submip against KP+feaspump on the same configs, since only the
   sub-problem solver differs — and it is what tells you whether the new heuristic is worth
   pursuing at all. The contender experiments are **explicitly not part of this work** (S8): they
   need the full published protocol (3600 s, 3 seeds, all three test sets, every competing method,
   Cochran's Q then pairwise McNemar with Holm correction) and they are a research exercise, not an
   implementation step.

   Practical consequence for whoever implements this: **finish at S7 and stop.** Report the baseline
   result, and do not treat the absence of S8 as unfinished work. If the baseline comes out flat,
   that is itself the useful finding and S8 may not be worth running.

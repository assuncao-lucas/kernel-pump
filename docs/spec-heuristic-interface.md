# Spec: a generic feasibility-heuristic interface for the Kernel Pump

**Status:** implemented and verified (M0-M5); M6 not done
**Branch:** `generic_interface` (from `main` @ `fc31aa97`)
**Scope:** `include/kernelpump/`, `src/kernelpump.cpp`, `src/main.cpp`, `tests/`, `settings/`

## 1. Goal

Today `KernelPump` owns a concrete `FeasibilityPump` by value and drives it bucket by bucket. The
method itself, however, is not specific to FP: KP needs *some* procedure that, given a sub-MIP
restricted to the active kernel, either returns an integer-feasible point or returns a fractional
point that is "closest to feasible" so KP can decide how to grow the kernel.

This spec defines an abstract **feasibility heuristic** interface, makes `FeasibilityPump` one
implementation behind it, and lets KP be configured to drive any other implementation — including
ones not yet written.

**Hard requirement:** with the default configuration (`kp.heuristic=feaspump`), every committed
config must produce the *same final behaviour* as today. Section 7 defines what "same" means
operationally and Section 9 says how it is verified.

## 2. Non-goals

- No change to the FP algorithm, nor to any `fp.*` key semantics.
- No change to kernel/bucket construction (`BuildKernelAndBuckets`), beyond removing one accidental
  dependency on FP (Section 6.3).
- No change to the `.sol` file format (see the constraint in Section 7.4).
- No second heuristic is *delivered* here. Appendix B sketches one only to prove the interface is
  sufficient; implementing it is follow-up work.
- Non-CPLEX backends stay as they are.

## 3. What KP actually requires of FP today

Extracted from every `feasibility_pump_` reference in `src/kernelpump.cpp` and
`include/kernelpump/kernelpump.h`. This is the complete surface — the interface needs nothing more.

| Call site | Member used | Purpose |
| --- | --- | --- |
| `kernelpump.cpp:818` | `readConfig()` | once, before the bucket loop |
| `kernelpump.cpp:892` | `init(model_)` → `bool` | per bucket; `false` ⇒ subproblem already infeasible |
| `kernelpump.cpp:910` | `pump(tl, stopNoImpr, xStartFrac, xStartDist, pFeas)` → `{found, progressed}` | per bucket; the actual work |
| `kernelpump.cpp:925` | `getSolution(v)` | on success |
| `kernelpump.cpp:926` | `getPrimalBound()` | on success |
| `kernelpump.cpp:927,951` | `getClosestDist()` | on success, and to test for a better reference point |
| `kernelpump.cpp:958` | `getClosestFrac(v)` | harvest the new best fractional point |
| `kernelpump.h:30` | `getIterations()` | reported as `# iterations` |
| `kernelpump.cpp:606,680` | `getReOptMethod()` | **not** heuristic work — KP's own kernel-repair LP solves (Section 6.3) |

Dead on arrival: `nbins()`, `bins()` and `getfirstOptMethod()` are declared in `feaspump.h` (with a
comment claiming `bins()` is "needed to build Kernel Pump's initial kernel and buckets") but have
**no caller anywhere** in `src/`, `include/`, `tests/` or `results/`. They do not belong in the
interface. Leave them on `FeasibilityPump` or delete them separately; do not generalise them.

## 4. The interface

New header `include/kernelpump/feasibility_heuristic.h`, in `namespace dominiqs`, following the
existing pluggable-component pattern (`SolutionTransformer` in `fp_interface.h`, `Ranker` in
`ranking.h`): abstract base + `SingletonHolder<Factory<...>>` + a static recorder.

```cpp
#ifndef FEASIBILITY_HEURISTIC_H
#define FEASIBILITY_HEURISTIC_H

#include <memory>
#include <string>
#include <vector>

#include <utils/factory.h>
#include <utils/singleton.h>

#include "kernelpump/mipmodel.h"

namespace dominiqs
{
	/** Warm start handed to a heuristic by its driver. */
	struct HeuristicStart
	{
		std::vector<double> frac;        /**< starting fractional point; empty ⇒ solve the relaxation */
		double dist = INFBOUND;          /**< driver's distance for `frac` (informational; see §7.3) */
		bool primalFeasible = false;     /**< is `frac` primal feasible? */

		bool empty() const { return frac.empty(); }
	};

	/** What one invocation of a heuristic achieved. */
	struct HeuristicOutcome
	{
		bool foundSolution = false;  /**< an integer-feasible point is available via getSolution() */
		bool progressed = false;     /**< the heuristic produced at least one reference point, i.e.
		                                  it was able to iterate at all. false ⇒ the driver should
		                                  treat the subproblem as (locally) infeasible. */
	};

	/**
	 * A procedure that looks for an integer-feasible point of the model it was init'd with, and
	 * that can report a "closest to feasible" fractional point when it fails.
	 *
	 * Lifecycle contract — see §5 of docs/spec-heuristic-interface.md. In short: ONE instance per
	 * driver run, readConfig() once, then init()+run() once per subproblem, with state carried
	 * across invocations at the implementation's discretion.
	 */
	class FeasibilityHeuristic
	{
	public:
		virtual ~FeasibilityHeuristic() = default;

		/** Read own parameters from gConfig(). Called once, before the first init(). */
		virtual void readConfig() {}

		/** Prepare for a subproblem. Returns false if the subproblem is already infeasible. */
		virtual bool init(MIPModelPtr model, const std::vector<char> &ctype = std::vector<char>()) = 0;

		/**
		 * Search for an integer-feasible point.
		 * @param timeLimit        wall-clock budget for this invocation
		 * @param start            warm start, or an empty HeuristicStart
		 * @param stopOnNoImprovement may the heuristic give up early on stagnation?
		 */
		virtual HeuristicOutcome run(double timeLimit, const HeuristicStart &start,
		                             bool stopOnNoImprovement) = 0;

		virtual bool foundSolution() const = 0;
		virtual void getSolution(std::vector<double> &x) const = 0;
		virtual double getPrimalBound() const = 0;

		/** Distance of the closest point to feasibility seen in the LAST invocation. */
		virtual double getClosestDist() const = 0;
		/** Projection of that closest point. Empty if none. */
		virtual void getClosestFrac(std::vector<double> &frac) const = 0;

		/** Total iterations over all invocations so far (reported as `# iterations`). */
		virtual int getIterations() const = 0;

		/** Name, for logging. */
		virtual std::string name() const = 0;
	};

	typedef std::unique_ptr<FeasibilityHeuristic> FeasibilityHeuristicPtr;
	typedef SingletonHolder<Factory<FeasibilityHeuristic, std::string>> HeuristicsFactory;

} // namespace dominiqs

#endif /* FEASIBILITY_HEURISTIC_H */
```

Notes on the shape:

- `pump()`'s five positional arguments (three of them defaulted) collapse into `HeuristicStart`.
  The tuple return becomes a named `HeuristicOutcome`; `std::get<1>(result)` at
  `kernelpump.cpp:916` is exactly the thing whose meaning nobody can guess from the call site.
- `run()` rather than `pump()`: "pump" names the FP mechanism, not the contract.
- No `resetPartial()` / `resetTotal()` in the interface. KP never calls either (`init()` calls
  `resetPartial()` internally); `main.cpp:331` calls `fp.resetTotal()` on the concrete FP in the
  standalone path. Keeping them off the interface prevents a future driver from assuming a reset
  semantics that implementations need not share.
- No LP-method accessor. See Section 6.3.

## 5. Lifecycle contract

This is the part that is easy to get wrong, because FP's behaviour across repeated invocations is
*stateful by design* and KP depends on that.

1. **One instance per KP run.** KP constructs the heuristic once (in `Init()`) and reuses it for
   every bucket. Constructing one per bucket is a behaviour change, not a refactor: `fp.iterLimit`,
   `fp.stage1IterLimit` and `fp.stage2IterLimit` are decremented by the iterations consumed at the
   end of every `pumpLoop` (`feaspump.cpp:1975-1985`) and are *not* restored by `init()` →
   `resetPartial()`. They are therefore **global budgets shared across all buckets**, and `nitr`
   accumulates across buckets. A fresh instance per bucket would hand every bucket a full budget
   and silently invalidate every published result.
2. **`readConfig()` exactly once**, before the first `init()`. For FP this is where the RNG is
   seeded (`rnd.setSeed(seed); rnd.warmUp()` at `feaspump.cpp:390-391`) and where the rounder is
   created and seeded. Calling it per bucket would re-seed the RNG each time.
3. **`init()` then `run()` per bucket, in that order.** `init()` may presolve afresh, so every
   vector crossing the boundary must be mapped (`presolveSolution` / `postsolveSolution`). If
   `init()` returns `false`, `run()` is not called.
4. **`getClosestDist()` / `getClosestFrac()` are per-invocation**, not cumulative: FP clears them in
   `resetPartial()` and resets `closestDist = INFBOUND` when given a warm start
   (`feaspump.cpp:691`). KP maintains its own running best (`closest_dist_`, `closest_frac_`) and
   feeds it back in. The interface keeps this split — do not "improve" it by making
   `getClosestDist()` cumulative.
5. **An implementation documents its own cross-invocation state.** The budget carry-over above is
   FP's choice, not a requirement of the interface. A new heuristic may be stateless across
   invocations; it must say so in its header.

## 6. Changes to KernelPump

### 6.1 Member and construction

```cpp
// include/kernelpump/kernelpump.h
#include "kernelpump/feasibility_heuristic.h"   // replaces #include "kernelpump/feaspump.h"
...
private:
    FeasibilityHeuristicPtr heuristic_;
    std::string heuristic_name_ = "feaspump";   // from kp.heuristic
    char kernel_repair_lp_method_ = 'S';        // from kp.kernelRepairLpMethod, see §6.3
```

`KernelPump::readConfig()` gains:

```cpp
heuristic_name_ = gConfig().get("kp.heuristic", std::string("feaspump"));
```

`KernelPump::Init()` creates it via the factory and fails loudly on an unknown name:

```cpp
heuristic_.reset(HeuristicsFactory::getInstance().create(heuristic_name_));
if (!heuristic_)
    throw std::runtime_error(fmt::format("Unknown kp.heuristic '{}'", heuristic_name_));
consoleLog("kpHeuristic: {}", heuristic_->name());
```

`KernelPump::getIterations()` becomes `heuristic_ ? heuristic_->getIterations() : 0` (it is
currently called from `main.cpp:258` unconditionally, including when `Init()` failed — today that
reads a default-constructed FP and yields 0; preserve that).

`KernelPump::Reset()` must **not** destroy `heuristic_`: `Reset()` is called from `Init()` before
the heuristic is created, and from `main.cpp:284` after the run. Leaving the pointer alone keeps
both call orders valid.

### 6.2 The bucket loop

A mechanical rewrite of `kernelpump.cpp:818-920`:

```cpp
heuristic_->readConfig();                                  // was feasibility_pump_.readConfig()
...
if (heuristic_->init(model_))
{
    bool stop_on_no_improvement = (curr_bucket_index != total_num_buckets - 1);

    HeuristicStart start;
    if (!reset_fp_initial_basis_at_new_loop_ && !closest_frac_.empty())
    {
        start.frac = closest_frac_;
        start.dist = closest_dist_;
        start.primalFeasible = true;
    }

    auto outcome = heuristic_->run(curr_time_limit_iteration, start, stop_on_no_improvement);
    found_int_feasible_solution = outcome.foundSolution;
    feasible_fp = outcome.progressed;
}
```

Everything downstream — the `first_bucket_to_iter_pump_` update, the success branch, and the five
kernel-absorption branches — stays **character-for-character as it is**. Rename the local
`feasible_fp` to `progressed` if you like, but change no condition, no ordering, and no comparison
(`lessThan(curr_fp_closet_dist, closest_dist_)` included).

### 6.3 Breaking the `getReOptMethod()` coupling

`BuildKernelAndBuckets` solves LPs of its own while repairing the initial kernel
(`kernelpump.cpp:606` and `:680`) and asks FP which simplex method to use. That is not heuristic
work, and a generic heuristic has no business answering it.

**There is a latent bug here that the refactor must not "fix" by accident.** `BuildKernelAndBuckets`
runs at `kernelpump.cpp:809`, but `feasibility_pump_.readConfig()` runs at `:818` — *after* it.
So `getReOptMethod()` at `:606`/`:680` always returns the constructor default
`DEF_REOPT_METHOD = 'S'` (`feaspump.cpp:149`), **never** the configured `fp.reOptMethod`.

Resolution: give KP its own key, defaulting to the literal `'S'`:

```cpp
// KernelPump::readConfig()
std::string m = gConfig().get("kp.kernelRepairLpMethod", std::string("default"));
// "default"->'S', "primal"->'P', "dual"->'D', "barrier"->'B', else throw
// (same mapping as feaspump.cpp:246-256)
```

and use `kernel_repair_lp_method_` at both call sites. Because every committed config leaves this
key unset, and the stale read always produced `'S'`, this is **behaviour-neutral for all 25
configs** while removing the coupling and making the choice explicit.

Do **not** instead wire these call sites to `fp.reOptMethod`, and do **not** move
`heuristic_->readConfig()` before `BuildKernelAndBuckets` while keeping an LP-method accessor on the
interface. Either would change the simplex method used during kernel repair for any config that
sets `fp.reOptMethod` to something other than `default` — a real behaviour change, invisible in the
diff.

## 7. Behaviour preservation

### 7.1 FeasibilityPump stays untouched

`FeasibilityPump` is adapted, not modified. Add `src/feaspump_heuristic.cpp` +
`include/kernelpump/feaspump_heuristic.h`:

```cpp
class FeasibilityPumpHeuristic final : public FeasibilityHeuristic
{
public:
    void readConfig() override { fp_.readConfig(); }
    bool init(MIPModelPtr model, const std::vector<char> &ctype = {}) override
        { return fp_.init(model, ctype); }
    HeuristicOutcome run(double timeLimit, const HeuristicStart &start, bool stopOnNoImprovement) override
    {
        auto r = fp_.pump(timeLimit, stopOnNoImprovement, start.frac, start.dist, start.primalFeasible);
        return {std::get<0>(r), std::get<1>(r)};
    }
    bool foundSolution() const override { return fp_.foundSolution(); }
    void getSolution(std::vector<double> &x) const override { fp_.getSolution(x); }
    double getPrimalBound() const override { return fp_.getPrimalBound(); }
    double getClosestDist() const override { return fp_.getClosestDist(); }
    void getClosestFrac(std::vector<double> &f) const override { fp_.getClosestFrac(f); }
    int getIterations() const override { return fp_.getIterations(); }
    std::string name() const override { return "feaspump"; }
private:
    FeasibilityPump fp_;
};
```

Registered with the established static-recorder idiom (cf. `transformers.cpp:289`,
`ranking.cpp:190`), which works because `libkp` is whole-archive linked:

```cpp
class HEUR_FACTORY_RECORDER
{
public:
    HEUR_FACTORY_RECORDER()
    {
#ifndef SILENT_EXEC
        std::cout << "Registering FeasibilityHeuristics...";
#endif
        HeuristicsFactory::getInstance().registerClass<FeasibilityPumpHeuristic>("feaspump");
#ifndef SILENT_EXEC
        std::cout << "done" << std::endl;
#endif
    }
};
HEUR_FACTORY_RECORDER my_heur_factory_recorder;
```

Why an adapter rather than `class FeasibilityPump : public FeasibilityHeuristic`: `git diff` on
`src/feaspump.cpp` and `include/kernelpump/feaspump.h` stays **empty**, so FP-path identity is
provable by inspection instead of only by testing. `feaspump.cpp` is 2481 lines, most of it
experimental branches that no paper config reaches; not editing it at all is worth one thin class.
The alternative is noted in Section 11.

### 7.2 Call-sequence identity

For `kp.heuristic=feaspump` the sequence of calls reaching `FeasibilityPump`, their arguments, and
their order relative to KP's own model mutations and LP solves must be identical to today. The
mock-based test in Section 9.2 asserts exactly this.

### 7.3 Two warm-start quirks to preserve, not repair

While adapting, you will notice both of these. Leave them.

- `pump()` sets `closestDist = INFBOUND` when handed a warm start (`feaspump.cpp:691`) and only
  *logs* `xStartDist`. So `start.dist` is informational — it does not seed FP's notion of the best
  distance. KP's `lessThan(curr_fp_closet_dist, closest_dist_)` test at `kernelpump.cpp:954` thus
  compares a from-scratch distance against KP's running best. That is current behaviour.
- with a warm start, `closestFrac` is assigned the start point immediately
  (`feaspump.cpp:686`), so `run()`'s `progressed` (`!closestFrac.empty()`, `feaspump.cpp:822`) is
  effectively always `true` once a warm start exists — the `!feasible_fp` absorption branch is
  reachable mainly on the first bucket or with `kp.resetFPBasisAtNewPump=1`.

Document both in `feasibility_heuristic.h` next to `HeuristicStart::dist` and
`HeuristicOutcome::progressed`, so a new implementation knows what the driver actually relies on.

### 7.4 The `.sol` format is positional — do not touch it

`results/src/main.cpp` parses solution files **by line offset**, not by label (e.g.
`for (int loop = 0; loop < 8; ++loop) getline(input, line);` at `results/src/main.cpp:286-288`,
and the same pattern at `:334`+). Inserting, removing or reordering any line in
`Solution::WriteToFile` silently shifts every field the aggregator reads.

Therefore: **no changes to `Solution` or `WriteToFile` in this work.** If a later heuristic needs a
new statistic, append it strictly after the last currently-written line and audit every positional
read in `results/` in the same change. The FP-specific field name `first bucket to iter pump` stays
as is; it is now understood as "first bucket where the heuristic managed to iterate".

## 8. Configuration

| Key | Default | Meaning |
| --- | --- | --- |
| `kp.heuristic` | `feaspump` | which registered heuristic KP drives |
| `kp.kernelRepairLpMethod` | `default` (→ `'S'`) | simplex method for the initial-kernel repair LPs (§6.3) |

Both unset in every existing config, so `settings/all/config*.cfg` need **no edits** and keep their
behaviour. An unknown `kp.heuristic` throws at `Init()` with the name in the message.

Each implementation owns its own config prefix and reads it in `readConfig()`; FP keeps `fp.*`
untouched. A new heuristic must not read `fp.*`.

`method=` in `main.cpp` is unchanged: `solver`, `feaspump`, `kernelpump`. The standalone
`method=feaspump` path keeps using a concrete `FeasibilityPump` (see M5 for the optional unification).

## 9. Test plan

The repo has no unit-test framework; `tests/grbmodel_test.cpp` establishes the pattern — a
standalone executable, hand-computed expectations, `check(cond, description)`, non-zero exit on
failure. Follow it.

### 9.0 First, measure determinism (prerequisite)

KP is **not** unconditionally deterministic: the per-bucket budget is
`time_left / (num_buckets + 1)` recomputed from the wall clock (`kernelpump.cpp:826`), and FP's
stopping criteria are partly time-based. Two runs of the same instance × config × seed on `main`
can legitimately differ.

So before writing any golden test, run the **unmodified** binary twice over the candidate set and
keep only the (instance, config, seed) triples whose `.sol` files agree on every non-timing field.
That set is the regression suite; the rest can only be compared statistically (9.3).

```bash
# on main, twice, into two folders, then diff ignoring the two timing lines
diff <(grep -vE '^(total time|time building kernel)' a/s_cfg_inst_1.sol) \
     <(grep -vE '^(total time|time building kernel)' b/s_cfg_inst_1.sol)
```

Pick a time limit generous enough that the method finishes well inside it — then the trajectory is
driven by the seeded RNG, not by the clock.

**Do not pick instances by size.** This was measured on `main @ fc31aa97` (CPLEX 2211, `config16`,
`seed=1`): of the 20 smallest MIPLIB instances in `instances/miplib.txt`, all 20 were perfectly
deterministic, but 17 of them reported `last bucket visited: 0` with 0–1 iterations — KP found a
solution on the initial kernel and **never entered the bucket loop**. They cannot detect a
regression in the code this work changes.

Select on *behaviour* instead, by probing candidates with a single run and keeping those that
satisfy both:

- `last bucket visited > 0` (the bucket loop actually ran), and
- `total time` comfortably below the time limit (not clock-bound).

A run that exhausts `fp.iterLimit` is ideal: it terminates early *and* deterministically, because
once the shared budget is gone the remaining buckets get zero iterations. `# iterations` equal to
`fp.iterLimit` exactly is the signature.

Probe results from `main` worth reusing as a starting pool (`config16`, `seed=1`, `tl=120`):

| instance | time | status | iters | buckets | last bkt |
| --- | --- | --- | --- | --- | --- |
| `supportcase18` | 0.5s | FOUND | 369 | 10 | 3 |
| `enlight_hard` | 0.5s | FAIL | 2005 | 11 | 2 |
| `tr12-30` | 0.07s | FOUND | 94 | 11 | 3 |
| `csched007` | 3.7s | FOUND | 4585 | 8 | 8 |
| `ns1830653` | 13.1s | FOUND | 1651 | 11 | 11 |
| `traininstance2` | 32.3s | FOUND | 11402 | 5 | 5 |
| `neos-3004026-krka` | 77s | FAIL | 12000 | 1 | 1 |

`peg-solitaire-a3` hit the 120s limit exactly and must be excluded — it is clock-bound by
construction.

### 9.1 Golden-file equivalence (the core requirement)

The golden set now lives in the repo as the **`kernel-pump-regression` skill**
(`.claude/skills/kernel-pump-regression/`): `verify.sh` replays every case and reports differences,
`golden/` holds the recorded `.sol` files and the manifest, `tools/` holds the capture/screen/probe
helpers. Note it does **not** live under `solutions/`, which holds committed research output that
the `results` binary rewrites in place (see the gotcha in the `kernel-pump-internals` skill), and
the 34 MB baseline binary is git-ignored rather than committed.

A flat "every config × every instance × 3 seeds" grid is needlessly expensive, because the slow
high-coverage instances are exactly the ones you least need repeated across seeds. Three tiers
instead:

| tier | content | grid | time limit |
| --- | --- | --- | --- |
| **A** | fast instances (< 2s), the regression workhorse | 4 KP configs × 14 instances × 3 seeds | 30s |
| **B** | slow instances with maximal bucket coverage | 4 KP configs × 4 instances × 1 seed | 120s |
| **C** | the paths M3 does not touch (3 FP + 2 solver configs) | 5 configs × 4 instances × 1 seed | 30s |

KP configs are `config16`, `config17`, `config21`, `config23`; FP are `config1`, `config3`,
`config20`; plain solver are `config18`, `config19` (the nine of
`settings/settings_list.txt`). Tier C exists only to prove the untouched paths stayed untouched, so
one seed suffices.

1. On `main` (or any commit where FP is still driven directly), capture all three tiers twice —
   run A is the golden set, run B is the determinism screen.
2. On `generic_interface`, rerun identically.
3. Compare field by field, ignoring `total time (s)` and `time building kernel and buckets (s)`.
   Everything else must match exactly: `STATUS`, `# iterations`, `# buckets`,
   `last bucket visited`, `first bucket to iter pump`, `value`, `reopt value`,
   `real integrality gap`, `projection integrality gap`, `num frac`, and all four `num bin vars *`
   fields.

`# iterations` and `last bucket visited` are the sharpest detectors here: the first catches any
accidental reset of FP's shared iteration budget (§5.1), the second catches any change in the
absorption branches.

Tier C covers the `method=feaspump` and `method=solver` paths, which M3 does not touch.

Reusable harness (written during M0, in the session scratchpad — move it somewhere durable if this
is repeated): `screen.sh` runs a grid twice and emits a per-case `deterministic YES/NO` verdict;
`capture.sh` captures one grid with a `MANIFEST` recording binary sha, commit, dirty files,
instance dir and time limit; `compare.sh` diffs two capture dirs field by field ignoring the two
clock fields and exits non-zero on any difference.

### 9.2 Interface-contract test (mock driver) — `tests/heuristic_interface_test.cpp`

Golden files prove the outcome; a mock proves the *mechanism*, including the branches no paper
config reaches. Register two test-only heuristics and run `KernelPump` on a small `.lp` written by
the test (CPLEX required, so gate the target on `CPLEX_FOUND`).

**`RecordingHeuristic`** wraps `FeasibilityPumpHeuristic`, forwards everything, and logs each call.
Assert:

- `readConfig()` called exactly once, before the first `init()`;
- `init()` and `run()` alternate, once per bucket, `#calls == last_bucket_visited + 1`;
- `run()`'s `stopOnNoImprovement` is `true` for every invocation except the last bucket's;
- `start.empty()` on the first invocation, and on later ones iff `kp.resetFPBasisAtNewPump=1`;
- when non-empty, `start.frac` equals KP's `closest_frac_` and `start.dist` its `closest_dist_`;
- `run()`'s `timeLimit` equals `time_left/(n+1)` for all but the last bucket, and the whole
  remaining budget for the last;
- exactly one heuristic instance is constructed per `KernelPump::Run` (count constructions);
- the same object is reused across buckets (compare `this`).

**`ScriptedHeuristic`** returns a canned `HeuristicOutcome` sequence with no solver work, so each
kernel-absorption branch can be driven deliberately. One case per row of the table in the
`kernel-pump-internals` skill:

| Scripted outcome | Expected kernel effect |
| --- | --- |
| `init()` → false | reference kernel absorbed; loop continues |
| `{false, false}` | whole reference kernel absorbed permanently |
| `{false, true}` + `forceBucketVarsIntoKernel=1` | whole reference kernel absorbed |
| `{false, true}` + improved closest frac + `force=0` | only binaries with `frac > 0` absorbed |
| `{false, true}` + no improvement + `force=0` | nothing absorbed; bucket vars leave |
| `{true, true}` | solution recorded, reference kernel absorbed, loop breaks |

Assert `getNumVarsInKernel()`, `getLastBucketVisited()` and `getFirstBucketToIterPump()` after each.
These branches are live code that **no committed config exercises** (all 25 set
`forceBucketVarsIntoKernel=1`), so this is the first coverage they get.

**Negative case:** `kp.heuristic=nope` throws from `Init()` with `nope` in the message.

### 9.3 Statistical check (only if 9.1 leaves gaps)

For instances that fail the 9.0 determinism screen, run the full batch on both branches and compare
success counts per config with `results/src/McNemar_test.py`. A significant difference is a red
flag. This is a weaker test — use it to cover what 9.1 cannot, not as a substitute.

### 9.5 Result of the verification (2026-11-20)

Baseline `fc31aa97`, CPLEX 2211, Release with `-DSILENT_EXEC`, golden set at
`.claude/skills/kernel-pump-regression/golden` (run it with that skill's `verify.sh`).

| run | outcome |
| --- | --- |
| M1+M2 only (adapter added, KP still driving FP directly) | 199/199 identical |
| M3+M4+M5 (KP driving the interface) | **191/191 identical, 0 differing** |
| `tests/heuristic_interface_test.cpp` | 45/45 checks pass |

Two process lessons, both of which produced false alarms before being understood:

1. **Two agreeing runs are not a determinism proof.** `neos-1456979` passed the 2-run M0 screen on
   several cases, then produced up to **3 distinct outcomes in 3 runs of the unmodified baseline**
   (config17: 3/3 distinct for every seed). All 8 of its cases were dropped, taking the set from
   199 to 191. Screen with at least 3 runs if the set is ever rebuilt.
2. **Never compile while verifying.** An earlier M3 run reported 2 differences; they came from
   rebuilding `build/kp` underneath the running comparison and from the CPU contention of a
   parallel build, not from the code. Re-run on an idle machine against a fixed binary: 0
   differences.

### 9.4 Build matrix

Configure and build on Linux/GCC with CPLEX, and with Gurobi both present and absent, to confirm the
new files are in `libkp` for every combination and that `results/` (which compiles
`../src/kernelpump.cpp`) still links.

## 10. Milestones

| # | Content | Done when |
| --- | --- | --- |
| **M0** | Determinism screen (9.0) + golden set captured on `main` | **DONE** — 204 cases run twice, 199 deterministic, archived with `verify.sh` |
| **M1** | `feasibility_heuristic.h`: interface, `HeuristicStart`, `HeuristicOutcome`, factory | **DONE** — header compiles; nothing else changed |
| **M2** | `FeasibilityPumpHeuristic` adapter + registration; `feaspump.*` untouched | **DONE** — empty `git diff` on `src/feaspump.cpp` and `include/kernelpump/feaspump.h`; factory resolves `"feaspump"` at runtime |
| **M3** | KP drives the interface (§6.1, §6.2) + `kp.heuristic` | **DONE** — code complete; golden-file equivalence (9.1) verified, see §9.5 |
| **M4** | `kp.kernelRepairLpMethod`, `getReOptMethod()` coupling removed (§6.3) | **DONE** — no `feasibility_pump_` identifier remains in `kernelpump.cpp` (landed with M3: the interface has no LP accessor, so M3 cannot compile without it) |
| **M5** | `tests/heuristic_interface_test.cpp` (9.2) | **DONE** — 45/45 checks pass |
| **M6** *(optional)* | `method=feaspump` in `main.cpp` also goes through the factory | `method=feaspump` golden files unchanged |

M0 genuinely comes first: without a golden set captured on `main`, M3 cannot be validated, and the
whole point of the exercise is that the FP path does not move.

## 11. Alternatives considered

- **`FeasibilityPump : public FeasibilityHeuristic` (no adapter).** Less indirection and the
  obvious OO move; the virtuals can be one-line inline forwarders in the header, so
  `feaspump.cpp` still needs no edit. Rejected as the primary route only because it perturbs
  `feaspump.h` and FP's vtable, and the adapter makes "FP is byte-identical" checkable with
  `git diff`. If the extra class later feels like dead weight, collapsing it is a mechanical
  follow-up.
- **A template parameter (`KernelPump<Heuristic>`) instead of a virtual interface.** Zero dispatch
  cost, but the choice must then be a compile-time constant, which defeats `kp.heuristic` as a
  config key. The dispatch happens once per bucket — tens of calls per run — so the virtual call is
  free in context.
- **Extending `pump()` with a strategy enum.** Keeps everything in `feaspump.cpp`, which is already
  2481 lines of mostly-unreachable experimental branches. That is the problem, not the solution.
- **A C-style callback** (cf. the unused `fp_c_interface.h`). No state across buckets, which §5.1
  requires.

## 12. Risks

| Risk | Mitigation |
| --- | --- |
| Accidental per-bucket reset of FP's shared iteration budget | §5.1; `# iterations` in 9.1 detects it; explicit mock assertion in 9.2 |
| `readConfig()` moved earlier, changing kernel-repair LP method | §6.3 removes the dependency entirely |
| Golden comparison too noisy to conclude anything | M0 screens for determinism before anything is compared |
| A new field added to `.sol` shifts `results`' positional parsing | §7.4: no `Solution` changes in this work |
| Factory registration stripped by the linker | same static-recorder idiom as transformers/rankers, which works because `libkp` is whole-archive linked |

---

## Appendix A: files touched

| File | Change |
| --- | --- |
| `include/kernelpump/feasibility_heuristic.h` | **new** — interface, structs, factory |
| `include/kernelpump/feaspump_heuristic.h` | **new** — FP adapter |
| `src/feaspump_heuristic.cpp` | **new** — adapter + factory registration |
| `include/kernelpump/kernelpump.h` | member becomes `FeasibilityHeuristicPtr`; include swap; two new params |
| `src/kernelpump.cpp` | `readConfig()`, `Init()`, bucket loop, two LP-method call sites |
| `CMakeLists.txt`, `results/CMakeLists.txt` | add `src/feaspump_heuristic.cpp` to `libkp`; add the test target under `CPLEX_FOUND` |
| `tests/heuristic_interface_test.cpp` | **new** — contract + branch tests |
| `src/feaspump.cpp`, `include/kernelpump/feaspump.h` | **none** (invariant of M2) |
| `src/solution.cpp`, `include/kernelpump/solution.h` | **none** (§7.4) |
| `settings/all/*.cfg` | **none** (§8) |
| `CLAUDE.md`, `.claude/skills/kernel-pump-internals/` | document `kp.heuristic` and the interface once M3 lands |

## Appendix B: sufficiency check — a second heuristic

Not part of this work; written only to confirm the interface needs nothing more. A **sub-MIP
heuristic** (`kp.heuristic=submip`) that mirrors `method=solver` on each kernel:

- `init(model)` — keep the model; return `model->presolve()`-style feasibility as FP does.
- `run(timeLimit, start, _)` — set `IntParam::SolutionLimit = 1`, optional MIP feasibility
  emphasis, `DblParam::TimeLimit = timeLimit`, inject `start.frac` as a MIP start if non-empty,
  `mipopt()`. On a solution: `{true, true}`. On proven infeasible: `{false, false}`. On
  time-out with no solution: `{false, true}`, with `getClosestFrac()` returning the LP relaxation
  point and `getClosestDist()` its integrality distance via `computeIntegralityGap`.
- `getIterations()` — node count, or 0.

Every member maps onto something `MIPModelI` already exposes, and the KP driver needs no change —
which is the actual test of this design. It would also make `kp.heuristic=submip` a meaningful
research baseline: kernel growth driven by an exact sub-solver instead of a pump.

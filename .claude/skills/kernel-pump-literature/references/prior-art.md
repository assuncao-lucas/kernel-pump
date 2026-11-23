# FP-family prior art

Per-paper notes on everything in `papers/`, written to answer one question at design time: **has
this been done, and what did it achieve?** Summaries only — read the PDFs for detail. All paths are
relative to the repository root.

Ordered roughly as the line of work developed.

---

## 1. The feasibility pump — the origin

`papers/feas-pump.pdf` · Fischetti, Glover & Lodi · *Math. Program., Ser. A* 104:91–104 (2005)

Defines FP. Keep two complementary trajectories: `x*` satisfies the linear constraints but may be
fractional; `x̃` is integer but may violate them. Each pumping cycle rounds `x*` to `x̃`, then solves
`min{Δ(x, x̃) : x ∈ P}` to get a new `x*`. Stop when `Δ = 0`. Random perturbation breaks cycling.

Tested on 83 hard 0-1 MIPs against CPLEX 8.1: competitive in speed and first-solution quality.
CPLEX found no feasible root-node solution on 19 instances; FP failed on only 3.

**For KP design:** this is the loop KP decomposes. Rounding is "fast and simple, but somehow blind"
— every later paper improves one of its steps.

---

## 2. Feasibility pump 2.0 — the direct ancestor of this code

`papers/feas-pump-2.0.pdf` · Fischetti & Salvagnin · *Math. Prog. Comp.* 1:201–222 (2009)

Replaces naive rounding with a **diving-like procedure based on rounding plus constraint
propagation** (a Constraint Programming tool). Substantially improves FP's success rate without
slowing it down, and markedly improves the quality of solutions found; the headline mechanism is a
**reduced number of iterations to convergence**.

**For KP design:** this is `fp.frac2int=propround`, which every paper config uses — the baseline
rounding, not an upgrade still available. `include/kernelpump/fp_c_interface.h` is literally
described as the "C interface to FP2.0" and is authored by Salvagnin, so this codebase descends from
that implementation. Proposals to "add constraint propagation to the rounding" are already in place.

---

## 3. Boosting the feasibility pump — line search in the rounding step

`papers/boost-feas-pump.pdf` · Boland, Eberhard, Engineer, Fischetti, Savelsbergh & Tsoukalas ·
*Math. Prog. Comp.* (2014)

Enhances the rounding step with a **clever integer line search** that efficiently explores a large
set of integer points, rather than committing to a single rounding of `x*`.

**For KP design:** "try more than one integer point per iteration" is taken. Note the difference in
kind: this enriches the *rounding*, KP decomposes the *variable set*. They are orthogonal and could
compose.

---

## 4. Multiple reference vectors and objective scaling

`papers/multi-vectors-feas-pump.pdf` · Mexi, Berthold & Salvagnin ·
*EURO J. Computational Optimization* 11:100066 (2023)

Observes that of 15+ FP papers, none used information **across** iterations — all optimised towards
a single reference integer point. Contributions:

- **mRENS**: use LP-feasible vectors from the main loop to tighten variable domains before the
  expensive enumeration stage (Stage 3).
- multiple integer reference vectors in Stages 1 and 2, via a conic combination of distance
  functions in the projection objective.
- alternative **objective scaling** terms balancing the distance functions against `cᵀx`.

Embedded setting: 29.3% solution-quality improvement, 4.0% faster, 16.0% fewer iterations. Code:
<https://github.com/GioniMexi/feaspump3> — and **FeasPumpCollection, which this repo extends, is the
same author's**.

**For KP design:** the most important entry here. Much of this is *already present but disabled* in
our `FeasibilityPump`: `fp.numIntegersObj`, `fp.numFracsObj`, `aggregateInts`, harmonic and
exponential-decay weights, `lessViolatedIntegers`, `bestObjIntegers`, the Stage-3/RENS family. So
"use multiple reference vectors" is neither novel nor unimplemented — the honest experiment is to
switch the existing keys on and measure, ideally combined with KP, which is also one of the KP
paper's own stated future directions (combine with OFP variants).

Their standalone protocol matches ours: 10 000 Stage-1 and 2 000 Stage-2 iterations. They tried 50
and 100 and report it "does not change the overall picture by much" — so iteration-limit tuning is a
known dead end.

---

## 5. Feasibility Jump — the LP-free competitor

`papers/feasibility-jump.pdf` · Luteberget & Sartor · *Math. Prog. Comp.* 15:365–388 (2023)

Not an FP variant: stochastic **guided local search over a Lagrangian relaxation**. Repeatedly pick
a variable and set it to minimise a weighted sum of constraint violations, raising the weights
(Lagrangian multipliers) of constraints still violated at local minima. Needs **no LP solution at
all**, so it can produce solutions before presolve finishes. Incomplete: no optimality, no
infeasibility detection.

Won 1st place at the MIP 2022 Computational Competition and now runs by default in FICO Xpress 9.0.
Deliberately not over-tuned, inviting further work; their own list of directions starts with richer
variable neighbourhoods.

**For KP design:** the relevant external baseline for "find a first feasible solution fast", and the
strongest argument that the LP-solve-per-iteration cost is FP's structural weakness. KP attacks that
cost by shrinking each LP (fewer active binaries); FJ attacks it by removing the LP. A KP variant
driven by an FJ-style inner heuristic is a natural experiment, and the `FeasibilityHeuristic`
interface makes it a drop-in rather than a rewrite.

---

## 6. Fix-propagate-repair — portfolio diving

`papers/fix-propagate-repair-matheuristic.pdf` · Salvagnin, Roberti & Fischetti ·
*Math. Prog. Comp.* (2024), DOI 10.1007/s12532-024-00269-5

Diving framework: fix-and-propagate **plus solution repair after each step**, the repair loosely
based on WalkSAT. Essentially LP-free — the full LP is solved only at the start (and only for
rankings that need it), with small LPs at the bottom of a dive for continuous variables. Different
variable-ranking and value-selection strategies give different heuristics.

The key finding: **individual strategies are not robust on a heterogeneous testbed, but a portfolio
is** — 189 of 240 MIPLIB 2017 benchmark instances in seconds.

**For KP design:** two things. First, `fp.walksatPerturbe` already exists here (disabled), so
WalkSAT-style repair is not new to this code. Second, their portfolio result is independent support
for the KP paper's own "parallel portfolio of KP variants" direction — and the KP⁻/KP⁺/KP\* variants
are already a natural portfolio, since they differ only in the kernel-driving objective.

---

## 7. Kernel search as a matheuristic

`papers/A kernel search heuristic for a fair facility location problem.pdf` · Filippi, Guastaroba,
Huerta-Muñoz & Speranza · *Computers & Operations Research* 132:105292 (2021)

Kernel Search applied to a fair facility location problem (minimising the conditional β-mean, between
p-median and p-center). Representative of the Kernel Search line: identify a promising "kernel" of
variables, sort the rest into buckets, and solve a sequence of restricted MIPs that grow the kernel.

**For KP design:** the other half of KP's ancestry. KP's contribution is transplanting this
decomposition **into FP** — so a proposed change to kernel or bucket construction should be checked
against the Kernel Search literature, where bucket sizing and ordering have been studied in their
own right, rather than reasoned about from scratch.

---

## Also in `papers/`, for instance provenance not method design

- `papers/traveling-tournament-problem.pdf` · Melo, Urrutia & Ribeiro · *J Scheduling* 12:607–622
  (2009) — the traveling tournament problem with predefined venues, the source of the **TTP-PV**
  test set (`instances/ttp-pv.txt`).
- `papers/artigo-joaquim.pdf` · Costa, Paniago, de Andrade, Noronha & Vieira · *Computer Networks*
  165:106939 (2019) — ILP formulations for variable data-rate/channel-bandwidth scheduling in
  wireless networks.

The **R-STOP-DP** set (`instances/r-stop-dp.txt`, robust Steiner team orienteering) is the other
structured test set; its defining paper is not in `papers/`.

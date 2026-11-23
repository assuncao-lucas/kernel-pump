---
name: kernel-pump-literature
description: The published record for this Feasibility Pump / Kernel Pump codebase — the Kernel Pump paper (Math. Prog. Comput. 2026) plus the FP-family prior art in papers/ — and what it implies for any proposed change to the method. Use BEFORE designing, proposing, specifying or implementing any algorithmic improvement to FP or KP, when judging whether an idea is novel or already tried, when choosing an experimental protocol or deciding what to compare against, when a config key's purpose is unclear, and when writing up results. Records which variants the authors already tested and rejected, which code paths are published versus experimental, and the exact protocol a new result must match to be comparable.
---

# The published record, and what it constrains

This codebase is the implementation behind a published paper. That has two consequences that are
easy to forget while reading C++:

1. **Some "odd" behaviour is the published method**, not a bug to fix. The shared iteration budget
   is the clearest example (see below).
2. **Some ideas that look new have already been tried and rejected by the authors**, and the config
   keys for them are still in the code. Re-proposing one is not an improvement.

Check a proposed change against this file *before* writing a spec or code.

## The paper of record

> **Kernel pump: Boosting the convergence of the feasibility pump through a kernel search-based
> decomposition.** Lucas Assunção, Sebastián Urrutia, Andréa Cynthia Santos.
> *Mathematical Programming Computation*, published 7 August 2026.
> <https://doi.org/10.1007/s12532-026-00333-2> — **open access, CC BY 4.0**.

Novelty claim, in the authors' words: *"To our knowledge, this is also the first study to introduce
a decomposition approach within the FP framework."* The core intuition is that **not all binary
variables are active (set to one) in a feasible solution**, so they can be decomposed kernel-search
style.

A June 2025 talk version exists at `~/Downloads/kernel_pump.pdf` under the older title *"Enhancing
the convergence of the Feasibility Pump through variable decomposition"* — useful for framing, but
the published text is authoritative and differs.

The FP-family prior art is in `papers/` (all tracked in git). Per-paper notes, including what each
one forecloses: `references/prior-art.md`.

## Variant names ↔ configs

The paper's variants map onto `settings/all/*.cfg` as follows (also in README.md). The superscript
encodes **which objective drives the kernel/initial LP**, not a different algorithm:

| Paper | Config | Objective used |
| --- | --- | --- |
| FP\* | `config1` | original MILP objective |
| FP⁻ | `config3` | minimise Σ binaries — the all-zeros extreme |
| FP⁺ | `config20` | maximise Σ binaries — the all-ones extreme |
| KP\* | `config16` | original MILP objective |
| KP⁻ | `config17` | minimise Σ binaries |
| KP⁺ | `config21` | maximise Σ binaries |
| KP⁺ᐟ⁻ | `config23` | mixed: kernel built maximising, FP's initial LP minimising |
| CPLEX_feas | `config18` | plain solve, feasibility emphasis |
| CPLEX_std | `config19` | plain solve, default |

## The published method, key by key

The paper's notation on the left, the code on the right. If a change alters any of these, it changes
the published method and must say so explicitly.

| Paper | Code |
| --- | --- |
| partition binaries into κ+1 ranges by distance to integrality `d(x̄ᵢ) = \|x̄ᵢ − round(x̄ᵢ)\|` | `kp.buildBucketsByRelaxationLayers=1` + `kp.sortByFractionalPart=1`, `kp.numBucketLayers=κ=10` |
| sub-divide each range by reduced-cost sign into Υ⁻ / Υ\* / Υ⁺; each non-empty sub-range is one bucket | the reduced-cost-sign split inside each layer (`src/kernelpump.cpp`) |
| kernel refinement: `min Σ_{i∉K} xᵢ` s.t. `x ∈ P`, absorbing every `x̃ᵢ > 0` | `kp.tryEnforceFeasibilityInitialKernel=1`, implemented with `CPXfeasopt` **MIN_SUM** relaxing the upper bounds of non-kernel binaries — the same objective by a different route, and the paper names "FeasOpt MinSum" |
| `K̃ ← K̃ ∪ Λ_λ` every iteration; no backtracking | `kp.forceBucketVarsIntoKernel=1` (monotone growth) |
| warm start each sub-problem with the least-fractional point so far, `frac(x̄)=Σ d(x̄ᵢ)` | `kp.resetFPBasisAtNewPump=0` |
| `δ ← δ_max/(\|Λ\|+1)`, remainder to the last bucket | hardcoded in `KernelPump::Run` |
| stage 1 ≤ 10 000 pumps, stage 2 ≤ 2 000 | `fp.stage1IterLimit=10000`, `fp.stage2IterLimit=2000` |
| ≤ 40 iterations without ≥10% improvement, reset per sub-problem **except the last bucket** | `fp.stage{1,2}NoImprIterLimit=40`, `fp.imprThreshold=0.9`, and `stopWithNoImprLimit` disabled for the final bucket |
| *"Maximum iteration count is shared globally and never reset"* | `fp.iterLimit=12000`, decremented across buckets and deliberately not restored |
| propagation-based rounding (from FP 2.0) | `fp.frac2int=propround` |
| pure distance objective, no stage 3 | `fp.alpha=0.0`, `fp.doStage3=0` |

**The shared global iteration budget is published, intended behaviour.** The
`kernel-pump-internals` skill calls it "the biggest trap" because it surprises readers of the code,
not because it is wrong. Do not "fix" it; changing it invalidates comparison with the paper.

## Already tried and rejected — do not re-propose

The authors report these as tested and not worth keeping. The keys still exist, so an idea can
accidentally reduce to one of them.

| Rejected idea | Key that turns it on | What the paper reports |
| --- | --- | --- |
| reverse the sign of the original objective when building the kernel | `kp.reverseObjectiveFunction=1` | *"did not perform well"*; reached feasibility only once it had degenerated back to plain FP |
| replace the objective with an empty (null) function | `kp.buildKernelBasedOnNullObjective=1` | same finding |
| order binaries by their raw LP relaxation value (non-decreasing or non-increasing) instead of by distance to integrality | `kp.sortByFractionalPart=0` with KP | same finding |
| discard the original objective altogether | the KP⁻/KP⁺ extremes, taken to their limit | success rate holds up but *"solution quality degrades substantially"* — a real trade-off, not a free win |

So: **an improvement built on a different kernel-ordering criterion needs an argument for why it is
not the already-rejected "sort by LP value"**, and anything that drops the objective must report
solution quality, not just success rate.

## Where KP is known not to help

- KP's benefit comes from problems where *"only a small proportion of the binary variables assume
  value one in a feasible solution."* Without that structure, the gain shrinks.
- On the heterogeneous MIPLIB 2017 set, KP's absolute success is higher (KP⁻ 204/233, KP⁺ 203/233 vs
  FP\* 197/233) but **not statistically significant**. The dramatic wins are on the structured sets:
  R-STOP-DP 48/48 vs FP\* 16/48, with KP⁻ averaging 12.44s against FP\*'s 3124.10s; TTP-PV similar,
  roughly halving time.
- Convergence metrics show high standard deviation — performance is inconsistent across problem
  types.

A change that only moves MIPLIB numbers is therefore hard to defend statistically; the structured
sets are where effects are detectable.

## Open directions the paper itself sanctions

From its future work:

- better initial-kernel and bucket construction;
- **parallel portfolio of different KP variants** (note the fix-propagate-repair paper found a
  portfolio far more robust than any single strategy — converging evidence);
- extending the decomposition idea to **other MILP heuristics** — which is exactly what the
  `FeasibilityHeuristic` interface now enables (`docs/spec-heuristic-interface.md`);
- combining KP with **objective-focused FP variants** (OFP).

## Experimental protocol a new result must match

To be comparable with the published tables:

| | |
| --- | --- |
| instances | MIPLIB 2017 (233, 7 infeasible excluded), R-STOP-DP (48), TTP-PV (36) — `instances/*.txt` |
| time limit | 3600 s |
| seeds | 3 per instance × algorithm |
| solver | CPLEX 22.11, **single thread** (`multiThreading=0`) |
| hardware | 64-bit AMD EPYC-Milan, 2.4 GHz, 100 GB RAM |
| statistics | Cochran's Q across algorithms, then pairwise McNemar with **Holm** correction at 5% — `results/src/McNemar_test.py` |

Report success rate *and* solution quality; the paper's own caveat about discarding the objective
exists because success rate alone is misleading.

## Is the idea novel? — check prior art first

`references/prior-art.md` has the detail. The short version of what is already taken:

| Idea | Already done by |
| --- | --- |
| better rounding via constraint propagation | FP 2.0 (Fischetti & Salvagnin, MPC 2009) — and it is what `propround` already is |
| integer line search inside the rounding step | Boosting FP (Boland et al., MPC 2014) |
| multiple reference vectors, objective scaling, mRENS bound tightening | Mexi, Berthold & Salvagnin (EJCO 2023) — **and partly implemented-but-disabled here** as `fp.numIntegersObj`, `fp.numFracsObj`, harmonic/exponential-decay weights |
| LP-free guided local search over a Lagrangian relaxation | Feasibility Jump (Luteberget & Sartor, MPC 2023) |
| fix-propagate-**repair** diving, WalkSAT-style repair, portfolio of strategies | Salvagnin, Roberti & Fischetti (MPC 2024) — note `fp.walksatPerturbe` already exists here, disabled |
| kernel search as a MIP matheuristic | Guastaroba/Speranza line; `papers/A kernel search heuristic for a fair facility location problem.pdf` |
| **kernel-search decomposition inside FP** | **this paper** — the contribution to protect |

Before claiming novelty, also check the vestigial FP options listed in `kernel-pump-internals`:
several "new ideas" are already coded here and merely switched off, which means the honest
experiment is to enable and measure them, not to reimplement them.

## Published vs experimental code paths

These exist in the code but are **not** part of the published method. Using one is fine; claiming
paper-comparable results while using one is not.

| Code path | Status |
| --- | --- |
| `kp.buildBucketsConsideringVariableDependency=1` | not in the paper (it has no dependency notion between buckets); no config in `settings/all` enables it |
| flat bucketing: `kp.buildBucketsByRelaxationLayers=0` with `kp.maxBucketSize` | not the published construction |
| `kp.forceBucketVarsIntoKernel=0` (selective absorption) | not published; the paper's loop grows the kernel unconditionally |
| `kp.heuristic`, `kp.kernelRepairLpMethod` | added after publication; defaults reproduce published behaviour exactly |
| FP stage-3 variants, analytic-centre FP, exponential/logistic distance objectives, WalkSAT perturbation | inherited from prior work, disabled in every paper config |

## Checklist before specifying a change

1. Is it in **Already tried and rejected**? If so, stop, or state what is different.
2. Is it already in **prior art**, or already coded and merely disabled? Then measure, don't invent.
3. Does it alter anything in **the published method, key by key**? Say so explicitly in the spec.
4. Does it depend on the structure KP exploits (few active binaries)? If not, expect small gains.
5. Which test set can actually show the effect? MIPLIB alone will likely not reach significance.
6. Will it change `.sol` output or published numbers? Then the regression suite
   (`kernel-pump-regression`) must pass, or the deviation must be deliberate and documented.
7. Does it touch solution quality as well as success rate? Report both.

## Related

- `kernel-pump-internals` — what the code actually does, and the invariants that break silently.
- `kernel-pump-regression` — proving a change did or did not alter behaviour.
- `docs/spec-heuristic-interface.md` — the pluggable-heuristic interface, which implements the
  paper's own "extend the decomposition to other MILP heuristics" direction.

# Config reference

Config files live in `settings/all/config*.cfg` and are listed in `settings/settings_list.txt`.
Format is one `key=value` per line. Any key can be overridden on the command line:

```bash
./build/kp instance.mps --config settings/all/config16.cfg seed=3 kp.maxBucketSize=20
```

Parsing is `mergeConfig(args, gConfig())` in `src/main.cpp:61`; command-line values win.

## Trap: inert keys are still present

**Every config file contains the full `kp.*` and `fp.*` key set, including keys that the selected
`method` ignores.** Reading `config1.cfg` you will see a complete kernel-pump block, but it is
`method=feaspump`, so none of it does anything. Likewise `config18`/`config19` are `method=solver`
and ignore *both* `kp.*` and `fp.*`.

Before concluding that a key matters for a run, check `method=` first.

| `method=` | Honours |
| --- | --- |
| `solver` | only `timeLimit`, `mipPresolve`, `mipFeasEmphasis`, `multiThreading`, `printSol` |
| `feaspump` | the above + all `fp.*`; **all `kp.*` inert** |
| `kernelpump` | the above + all `fp.*` + all `kp.*` |

## Top-level keys

| Key | Default | Meaning |
| --- | --- | --- |
| `method` | *(none — throws)* | `solver` \| `feaspump` \| `kernelpump`. Case-insensitive. |
| `solver` | `cpx` | Backend. Only `cpx` works; others throw unless compiled in. |
| `runName` | `default` | Goes into the output filename `s_<runName>_<prob>_<seed>.sol`. By convention the config name. |
| `seed` | `0` | Seeds the rounder and the ranker. The only source of run-to-run variation. |
| `timeLimit` | `1e20` | Global wall-clock budget in seconds, split across buckets by KP. |
| `solutionFolder` | `../solutions/test/` | Output directory. Must already exist. |
| `mipPresolve` | `1` | Presolve before pumping. See the presolve-mapping section of the main skill. |
| `mipFeasEmphasis` | `0` | `method=solver` only: CPLEX feasibility emphasis. |
| `multiThreading` | `0` | When 0, forces `Threads=1`. Experiments run single-threaded. |
| `modelLogging` | `0` | Solver log. Suppressed entirely under `-DSILENT_EXEC`. |
| `printSol` | `0` | Dump nonzero variable values to stdout. |
| `testset` | `unknown` | Free-form label, logged only. |

## `kp.*` — Kernel Pump

Read in `KernelPump::readConfig` (`src/kernelpump.cpp`).

### Kernel/bucket construction

| Key | Default | Meaning |
| --- | --- | --- |
| `kp.maxBucketSize` | `100` | Flat mode: size of the kernel and of each bucket. Ignored in layered mode. |
| `kp.numBucketLayers` | `10` | Layered mode: number of equal-width value layers. Realised bucket count is usually larger (layers subdivide by reduced-cost sign) and empty groups are dropped. |
| `kp.buildBucketsByRelaxationLayers` | `0` | `1` = layered mode, `0` = flat mode. |
| `kp.buildBucketsConsideringVariableDependency` | `0` | Also pull in binaries co-occurring in a constraint with LP value > 0. Costly: builds an O(n²) co-occurrence bitset matrix. |
| `kp.sortByFractionalPart` | `0` | `1` = order binaries by `\|round(v) − v\|` (most fractional last); `0` = order by raw LP value, descending. |
| `kp.tryEnforceFeasibilityInitialKernel` | `0` | Iteratively enlarge the initial kernel via CPLEX `feasopt` until its LP is feasible. |

### Objective used to *build* the kernel

These replace the objective of the cloned LP whose solution defines the ordering. Mutually
exclusive in practice; if none is set, the instance's own objective is used.

| Key | Default | Meaning |
| --- | --- | --- |
| `kp.buildKernelBasedOnNullObjective` | `0` | Objective ≡ 0 (pure feasibility; ordering then driven by reduced costs). |
| `kp.buildKernelBasedOnSumVarsObjective` | `0` | Objective = Σ binaries. |
| `kp.buildKernelBasedOnSumVarsObjectiveMaxSense` | `0` | Modifier on the above: MAX instead of MIN. |
| `kp.reverseObjectiveFunction` | `0` | Flip whatever sense resulted. |

### Kernel growth

| Key | Default | Meaning |
| --- | --- | --- |
| `kp.forceBucketVarsIntoKernel` | `0` | `1` = absorb the whole bucket into the kernel every iteration (monotone growth). `0` = absorb only binaries with `closest_frac > 0`. **All paper configs use 1.** |
| `kp.resetFPBasisAtNewPump` | `0` | `1` = each bucket's FP solves its own initial LP. `0` = warm-start from the best fractional point so far. |

### Which heuristic KP drives

| Key | Default | Meaning |
| --- | --- | --- |
| `kp.heuristic` | `feaspump` | Name of the `FeasibilityHeuristic` implementation KP drives per bucket, resolved through `HeuristicsFactory`. `feaspump` is the classical FP via the `FeasibilityPumpHeuristic` adapter. An unregistered name is rejected early in `src/main.cpp` (loudly, before the try/catch that `-DSILENT_EXEC` would mute) and the message lists the registered names. No config in `settings/all` sets this key. |
| `kp.kernelRepairLpMethod` | `default` | Simplex method for the LPs solved while repairing the initial kernel (`kp.tryEnforceFeasibilityInitialKernel=1`). Same spellings as `fp.reOptMethod`: `default`→`'S'`, `primal`→`'P'`, `dual`→`'D'`, `barrier`→`'B'`. **Historical note:** these LPs used to call `FeasibilityPump::getReOptMethod()` *before* FP's `readConfig()` had run, so they always got the constructor default `'S'` regardless of `fp.reOptMethod`. This key makes the choice explicit; its default reproduces that `'S'` exactly, so no committed config changes behaviour. |

## `fp.*` — Feasibility Pump

Read in `FeasibilityPump::readConfig` (`src/feaspump.cpp:229`).

### Core

| Key | Default | Meaning |
| --- | --- | --- |
| `fp.frac2int` | `propround` | Rounder: `propround` (round-one-then-propagate) or `std` (independent threshold rounding). |
| `fp.ranker` | `FRAC` | Order in which `propround` rounds: `FRAC` (least fractional first, binaries before general integers), `LR`, `RND`. |
| `fp.alpha` | `0.0` | Initial weight of the original objective vs. the distance. `0.0` = pure distance (classical FP). |
| `fp.alphaFactor` | `0.9` | Per-iteration decay of alpha; snaps to 0 below `1e-4`. |
| `fp.integralityEps` | `1e-6` | Integrality tolerance. Overwritten at startup by the solver's own tolerance (`src/main.cpp:153`). |
| `fp.filterConstraints` | `1` | Drop numerically nasty rows from the propagator (all-continuous with coefficient dynamism > 10, or any row with dynamism > 1000). |
| `fp.randomizedRounding` | `1` | Randomise the rounding threshold instead of using 0.5. |

### Budgets and stopping

| Key | Default | Meaning |
| --- | --- | --- |
| `fp.iterLimit` | `100` | **Global** iteration budget. Under KP it is decremented across buckets, not reset — see the budget carry-over section of the main skill. |
| `fp.stage1IterLimit` | `50` | Stage-1 cap; also decremented across KP buckets. |
| `fp.stage2IterLimit` | `50` | Stage-2 cap; likewise. |
| `fp.stage1NoImprIterLimit` | `70` | Abort stage 1 after this many iterations without sufficient improvement. |
| `fp.stage2NoImprIterLimit` | `600` | Same for stage 2. |
| `fp.imprThreshold` | `0.9` | "Sufficient improvement" means new distance < `imprThreshold ×` best, so `0.9` ⇒ 10%. |
| `fp.doStage3` | `0` | Enable the stage-3 sub-MIP. Off in all paper configs. |
| `fp.lpIterMult` | `-1.0` | If > 0, caps simplex iterations per LP at `lpIterMult ×` root LP iterations (floor 10). Disabled by default. |

### Initial-LP objective (mirrors the `kp.*` kernel-objective block)

| Key | Default | Meaning |
| --- | --- | --- |
| `fp.forceNullObjectiveInitialLP` | `0` | Objective ≡ 0 for the initial LP. |
| `fp.forceSumVarsObjectiveInitialLP` | `0` | Objective = Σ binaries for the initial LP. |
| `fp.forceSumVarsObjectiveInitialLPMaxSense` | `0` | Modifier: MAX instead of MIN. |
| `fp.reverseObjectiveFunction` | `0` | Flip the resulting sense. |

### LP method

`fp.firstOptMethod` / `fp.reOptMethod` ∈ {`default`→`S`(simplex), `primal`→`P`, `dual`→`D`,
`barrier`→`B`}; an unknown value throws. Note KP's kernel-building LP always uses **dual**
simplex regardless of these (hardcoded, `src/kernelpump.cpp:264`), because barrier often returns
an optimal value without a usable primal solution.

### Experimental (not used by any paper config)

Present in the header and `readConfig`, but off in every config under `settings/all` that the
article uses. Confirm reachability before investing time in these branches:

`numIntegersObj`, `numFracsObj` (>1 ⇒ multiple reference vectors), `aggregateInts`,
`fracScaleFactor`, `intScaleFactor`, `harmonicWeights`, `exponDecayWeights`,
`lessViolatedIntegers`, `lessDistanceIntegers`, `bestObjIntegers`, `expObj`, `logisObj`,
`analcenterFP`, `walksatPerturbe`, `randomizeLP`, `penaltyObj`, `newScaleC`, `newScaleDelta`,
`avgFlips`, and the stage-3 family (`normalMIPStage3`, `newStage3`, `rensStage3`,
`rensClosestDistStage3`, `multirensStage3`, `stage3*`). Also `pdlpTol`,
`pdlpTolDecreaseFactor`, `pdlpWarmStart` — PDLP backend only, which is not production-ready.

## Config number → article name

From `README.md`, cross-checked against the config files.

| Config | Article | `method` | Kernel objective | FP initial-LP objective |
| --- | --- | --- | --- | --- |
| `config1` | FP\* | feaspump | — | original |
| `config3` | FP⁻ | feaspump | — | MIN Σx |
| `config20` | FP⁺ | feaspump | — | MAX Σx |
| `config16` | KP\* | kernelpump | original | original |
| `config17` | KP⁻ | kernelpump | MIN Σx | MIN Σx |
| `config21` | KP⁺ | kernelpump | MAX Σx | MAX Σx |
| `config23` | KP⁺ᐟ⁻ | kernelpump | MAX Σx | MIN Σx |
| `config24` | KP⁰ | kernelpump | null | null |
| `config19` | CPLEX_std | solver | — | — |
| `config18` | CPLEX_feas | solver (`mipFeasEmphasis=1`) | — | — |

So the superscript decodes as the objective driving kernel construction / the initial LP:
`*` = the instance's own objective, `⁻` = minimise Σ binaries, `⁺` = maximise Σ binaries,
`⁺ᐟ⁻` = maximise for the kernel but minimise inside the FP, `⁰` = null.

All KP configs above additionally share: `sortByFractionalPart=1`,
`tryEnforceFeasibilityInitialKernel=1`, `buildBucketsByRelaxationLayers=1`, `numBucketLayers=10`,
`maxBucketSize=10`, `forceBucketVarsIntoKernel=1`, `resetFPBasisAtNewPump=0`,
`buildBucketsConsideringVariableDependency=0`.

All FP and KP configs share: `fp.frac2int=propround`, `fp.alpha=0.0`, `fp.iterLimit=12000`,
`fp.stage1NoImprIterLimit=40`, `fp.stage2NoImprIterLimit=40`, `fp.imprThreshold=0.9`,
`fp.doStage3=0`, `fp.numFracsObj=1`, `fp.firstOptMethod=default`, `fp.reOptMethod=default`,
`timeLimit=3600`, `mipPresolve=1`.

`TESTSET.md` documents the remaining numbered configs (config2, 4–15, 22, 25), which cover
intermediate variants and the multiple-reference-vector experiments. `config25` is `config16`
with the objective reversed. Note `settings/all/config12cfg` is missing its `.` — it will not
match a `*.cfg` glob.

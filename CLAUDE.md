# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A C++17 research codebase implementing the **Feasibility Pump (FP)** family of MIP heuristics and an improved
variant called the **Kernel Pump (KP)**, built by extending the public FeasPumpCollection. It is used to run
large batches of experiments (instance × config × seed) over MIPLIB-style instances and produce the
result tables/statistics for an accompanying research article.

## Build

### Requirements

| Requirement | Notes |
| --- | --- |
| C++17 compiler | GCC 15 works. Headers must include `<cstdint>`/`<string>` explicitly — newer libstdc++ no longer pulls them in transitively. |
| CMake ≥ 3.10 | CMake 4.x removed compatibility with `cmake_minimum_required` below 3.5, so all subprojects (including vendored `extern/fmt`) declare ≥ 3.10. |
| Boost headers (`libboost-dev`) | **Header-only**, no compiled Boost libs are linked. Only `boost/dynamic_bitset.hpp` and `boost/iterator/iterator_facade.hpp` are used. A missing `libboost-dev` is the most common cause of a failing fresh clone. |
| Eigen3 | Located via `find_package(Eigen3)`, falling back to the `EIGEN3_INCLUDE_DIR` cache variable (defaults to `/usr/include/eigen3` on Unix). Override with `-DEIGEN3_INCLUDE_DIR=...`. |
| CPLEX | The only functional solver backend. Pass `-DCPLEX_ROOT_DIR=...`. |

On Debian/Ubuntu: `sudo apt install libboost-dev libeigen3-dev cmake build-essential`.

Two independent CMake projects share the same solver-detection machinery: the root project (`kp` executable +
`libkp` static lib) and `results/` (a separate `results` executable for aggregating solution files into tables).

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -S=.. \
  -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio2211/cplex \
  -DXPRESSDIR=/opt/fico/xpressmp87 \
  -DSCIP_DIR=/opt/scip
make -j12
cd ..
```

**Only the CPLEX interface (`solver=cpx`) is currently functional**; XPRESS/SCIP/PDLP model classes exist
(`src/xprsmodel.cpp`, `src/scipmodel.cpp`, `src/pdlpmodel.cpp`) and compile conditionally (`HAS_CPLEX`,
`HAS_XPRESS`, `HAS_SCIP`, `HAS_ORTOOLS` feature macros set by CMake `find_package` results), but are not
production-ready. Solver availability is auto-detected — omit a `-D..._DIR` if that solver isn't installed and
its code paths just won't compile in.

By default the build defines `-DSILENT_EXEC` (set in the top of `CMakeLists.txt`), which suppresses solver
logging for fast batch runs. Comment out that `add_definitions(-DSILENT_EXEC)` line and rebuild for verbose
solver output when debugging.

The results-aggregation subproject builds the same way from its own directory:

```bash
cd results && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -S=.. -DCPLEX_ROOT_DIR=... -DXPRESSDIR=... -DSCIP_DIR=...
make -j12
cd ../..
```

Eigen3 is resolved via `find_package(Eigen3)` with an `EIGEN3_INCLUDE_DIR` cache-variable fallback (see the
bottom of both `CMakeLists.txt` files). `libkp`, `results`, and `kp` are all whole-archive-linked against the vendored
`extern/fmt`, `extern/utils` (`Utils::Lib`), and `extern/propagator` (`Prop::Lib`) subprojects, which are pulled
in via `add_subdirectory` only when this is the master (top-level) CMake project.

Git branch/commit are baked into the binary at configure time via `include/kernelpump/version.h.in` →
`version.h` (`KP_GIT_HASH`, `KP_GIT_BRANCH`), so re-run `cmake` (not just `make`) after committing if you need
the embedded hash to be current.

## Running

```bash
./build/kp PATH_INSTANCE_FILE --config CONFIG_FILE seed=INTEGER solutionFolder=PATH timeLimit=FLOAT
```

Any config key can be overridden inline: `./build/kp instance.mps -c config.cfg kp.maxBucketSize=20`.

`method=` (in the config or as an override) selects the code path in `src/main.cpp`: `solver` (plain
CPLEX MIP solve), `feaspump`, or `kernelpump`. Config files live in `settings/all/*.cfg` and are grouped/listed
in `settings/settings_list.txt`. `fp.*` keys configure `FeasibilityPump`, `kp.*` keys configure `KernelPump`;
see any file under `settings/all` for the full key set. `TESTSET.md` documents what each numbered config
(config1…config25) varies, and the README maps the article's FP*/FP-/FP+/KP*/KP-/KP+ naming onto specific
config numbers.

There is no unit test suite. Correctness is checked at runtime: after `kp` finds a solution it re-verifies
every constraint against it (see the feasibility double-check block near the end of `src/main.cpp`) and throws
if any constraint is violated. Statistical validation across configs (e.g. McNemar's test) is done separately
in `results/src/McNemar_test.py` against `results/build`-generated CSVs — requires `statsmodels` and `pandas`.

## Running full experiment batches

Three shell scripts drive batch experiments over instance × config × seed grids; they all `cd build && make`
before running so the binary is always current:

- `script_run_experiments.sh --instances-dir PATH --instances-list FILE --solutions-dir PATH --configs-dir PATH --configs-list FILE --time-limit FLOAT --num-seeds INT [--percentage-random-sample INT]` — the general-purpose batch runner. Supports randomly sampling a percentage of the instance list.
- `script_reproduce_experiments.sh --instances-dir PATH --instances-list FILE --configs-list FILE --time-limit FLOAT --num-seeds INT [--percentage-random-sample INT]` — reproduces the paper's experiments; auto-compiles in silent mode and writes results into a timestamped subfolder of `solutions/`.
- `script_generate_results.sh` — post-processes a solutions folder into `summary.csv` and the LaTeX table bodies used in the article (via the `results` executable).
- `script` (older/ad-hoc variant of the runner, has instance/solver paths hardcoded near the top — treat as a template, not a stable entry point).

Instance lists referenced by these scripts live in `instances/*.txt` (e.g. `instances/miplib.txt`,
`instances/r-stop-dp.txt`, `instances/ttp-pv.txt`) — one relative instance path per line. Output solution files
land in `solutions/<run>/s_<config>_<instance>_<seed>.sol`, alongside `sorted_files_list.txt` (the resolved
instance list actually used, useful for reproducing random samples) and `summary.csv` (one row per
instance × config × seed with `status` and `total time` as the key columns).

## Architecture

`src/main.cpp` is the single CLI entry point. It parses args/config (`gConfig()`, `dominiqs` utils), constructs
a solver-specific `MIPModelPtr` (`CPXModel`/`XPRSModel`/`SCIPModel`/`PDLPModel`, all implementing the abstract
`MIPModelI` interface in `include/kernelpump/mipmodel.h`), reads the instance, and then dispatches to exactly
one of three code paths based on `method=`:

- **`solver`** — plain MIP solve via the underlying solver, stopped at the first feasible solution.
- **`feaspump`** — `FeasibilityPump` (`src/feaspump.cpp` / `include/kernelpump/feaspump.h`): the classical FP
  loop (LP relax → round → re-optimize) with pluggable rounding via `transformers.cpp` and variable-selection
  strategies via `ranking.cpp`.
- **`kernelpump`** — `KernelPump` (`src/kernelpump.cpp` / `include/kernelpump/kernelpump.h`): builds an initial
  "kernel" of variables plus a sequence of "buckets" of remaining binary variables (grouped by relaxation
  layers or variable dependency — see `buildBucketsByRelaxationLayers`/`buildBucketsConsideringVariableDependency`
  in `TESTSET.md`/config files), then incrementally grows the kernel bucket-by-bucket, delegating to an
  internal `FeasibilityPump` instance (`feasibility_pump_`) at each step to find/repair integer feasibility.

After a solution is found, `main.cpp` fixes the integer variables and re-solves the LP relaxation to get a
"reoptimized" objective value, then writes results via `Solution::WriteToFile` (`src/solution.cpp` /
`include/kernelpump/solution.h`) to `<solutionFolder>/s_<runName>_<probName>_<seed>.sol`.

The `results/` subproject is a standalone consumer: `results/src/main.cpp` + `results/src/csv_reader.hpp` walk
a solutions folder and `results/bestKnownBounds.csv` to compute gaps/summaries and emit the LaTeX table bodies
under `results/tables/`. `results/src/McNemar_test.py` runs pairwise statistical significance tests (McNemar,
with multiple-testing correction) across configs from that output.

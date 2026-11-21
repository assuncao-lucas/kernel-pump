# Kernel Pump and Feasibility Pump Collection

At this moment, the code supports the Feasibility Pump (FP) versions [[1]](#1),[[2]](#2),[[3]](#3),[[4]](#4),[[5]](#5) and the improved method called Kernel Pump (KP) [[6]](#6).
This code was developed by adapting and extending the freely available [FeasPumpCollection](https://github.com/GioniMexi/FeasPumpCollection).

## Publication

This repository is the implementation accompanying the following open-access article, which
describes the Kernel Pump method and reports the computational experiments reproducible with the
scripts below:

> L. Assunção, S. Urrutia, and A. C. Santos. **Kernel pump: Boosting the convergence of the
> feasibility pump through a kernel search-based decomposition.** *Mathematical Programming
> Computation*, 2026. [doi:10.1007/s12532-026-00333-2](https://doi.org/10.1007/s12532-026-00333-2)

Published under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

If you use this code, please cite:

```bibtex
@article{AssuncaoUrrutiaSantos2026,
  author  = {Assun{\c{c}}{\~a}o, Lucas and Urrutia, Sebasti{\'a}n and Santos, Andr{\'e}a Cynthia},
  title   = {Kernel pump: Boosting the convergence of the feasibility pump through a kernel
             search-based decomposition},
  journal = {Mathematical Programming Computation},
  year    = {2026},
  doi     = {10.1007/s12532-026-00333-2}
}
```

## Code overview

All the code is available at the folder [src](src), being the file [src/main.cpp](src/main.cpp) the starting point. The main KP and FP codes are in [src/kernelpump.cpp](src/kernelpump.cpp) and [src/feaspump.cpp](src/feaspump.cpp). Interfaces to the supported LP solvers are in [src/cpxmodel.cpp](src/cpxmodel.cpp), [src/xprsmodel.cpp](src/xprsmodel.cpp), [src/scipmodel.cpp](src/scipmodel.cpp), and [src/pdlpmodel.cpp](src/pdlpmodel.cpp).
Transformers (objects responsible for rounding a solution) are in [src/transformers.cpp](src/transformers.cpp). Different strategies for sorting variables before rounding the FP are in [src/ranking.cpp](src/ranking.cpp).

## Compilation instructions

This project uses CMake. For compiling, perform the following steps on the root folder of the project. Use the correct paths on your machine for CPLEX, XPRESS, SCIP:

- mkdir build
- cd build
- cmake -DCMAKE_BUILD_TYPE=Release -S=.. -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio2211/cplex -DXPRESSDIR=/opt/fico/xpressmp87 -DSCIP_DIR=/opt/scip..
- make -j12
- IF YOU NEED TO RUN IN **VERBOSE MODE**, comment line 'add_definitions(-DSILENT_EXEC)' in the main [CMAKELists.txt](CMakeLists.txt) file **BEFORE** executing 'make -j12'.
- cd ..

For now, only the CPLEX interface is ready for use.

Dependencies: a C++17 compiler, CMake ≥ 3.10, Boost headers and zlib (on Debian/Ubuntu: `sudo apt install build-essential cmake libboost-dev zlib1g-dev`).

### Windows

Requirements: [Visual Studio Build Tools 2022](https://visualstudio.microsoft.com/downloads/) (workload "Desktop development with C++"), CMake, Ninja and [vcpkg](https://github.com/microsoft/vcpkg) with the environment variable `VCPKG_ROOT` pointing to it. Boost and zlib are declared in [vcpkg.json](vcpkg.json) and installed automatically. From a "Developer PowerShell for VS 2022":

- cmake --preset windows-release -DCPLEX_ROOT_DIR="C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex"
- cmake --build build

The results sub-project builds the same way: `cmake -S results -B results/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=%VCPKG_ROOT%/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows-static-md -DCPLEX_ROOT_DIR=...` followed by `cmake --build results/build`.

In VS Code, install the C/C++ and CMake Tools extensions and pick the `windows-release` preset. The experiment scripts (`*.sh`) run under Git Bash.

### Gurobi (experimental)

If `GUROBI_HOME` is set (or `-DGUROBI_HOME=...` is passed), CMake also builds the Gurobi backend. Select it with `solver=grb` and always name the license file explicitly, e.g.:
```
./kp instance.mps -c settings/all/config16.cfg solver=grb grb.licenseFile=/path/to/gurobi.lic
```
Without `grb.licenseFile` the run stops before starting Gurobi. MIP presolve is not yet supported with Gurobi, so FP/KP run on the original model.

### Results

There is a sub-project dedicated to generating the results summary and LaTeX tables. It should also be compiled similarly as the main project. 

- cd ./results
- mkdir build
- cd build
- cmake -DCMAKE_BUILD_TYPE=Release -S=.. -DCPLEX_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio2211/cplex -DXPRESSDIR=/opt/fico/xpressmp87 -DSCIP_DIR=/opt/scip..
- make -j12
- cd ../..

## Usage

The main execution file is stored in the folder [build](build) after compiled.
```
./kp PATH_INSTANCE_FILE --config (-c) CONFIG_FILE seed=INTEGER solutionFolder=PATH timeLimit=FLOAT
```
Seed can be any integer number, and all the directories and files given as argument should exist.

Is is also possible to overrride config arguments like this:
```
.kp PATH_INSTANCE_FILE --config (-c) CONFIG_FILE arg=value
```

Check examples of config files in the folder [settings](settings/all).

## Experiments reproduction

In the main experiments of the Kernel Pump method, three benchmark of instances were used:
- [MIPLIB 2017 benchmark](https://miplib.zib.de/downloads/benchmark.zip)
- [Robust Steiner Team Orienteering Problem instances](https://drive.google.com/file/d/1LPPfpt_mbNgHhu0PWz8pNHT4-Qvi9ACa/view?usp=drive_link)
- [Traveling Tournament Problem with Predefined Venues](https://drive.google.com/file/d/1cf_0n4XnAH7WJKVh6aOf6OI10p2Ep7c0/view?usp=drive_link)

The specific instances used from each benchmark are listed in the files:
- [instances/miplib.txt](instances/miplib.txt)
- [instances/r-stop-dp.txt](instances/r-stop-dp.txt)
- [instances/ttp-pv.txt](instances/ttp-pv.txt)

The configurations tested are stored in the [settings/all](settings/all) folder, and listed in [settings/setting_list.txt](settings/settings_list.txt). The raw name of each configuration tested in the article is displayed as follows:
- config1 corresponds to FP<sup>*</sup>
- config3 corresponds to FP<sup>-</sup>
- config20 corresponds to FP<sup>+</sup>
- config16 corresponds to KP<sup>*</sup>
- config17 corresponds to KP<sup>-</sup>
- config21 corresponds to KP<sup>+</sup>
- config23 corresponds to KP<sup>+/-</sup>
- config19 corresponds to CPLEX<sub>std</sub>
- config18 corresponds to CPLEX<sub>feas</sub>

It is possible to run all the experiments by calling the script [script_reproduce_experiments.sh](script_reproduce_experiments.sh). In this case, **the code is automatically compiled to run in silent mode**, as to avoid unnecessary printing time in the final executions.

Script Usage:
```
./script_reproduce_experiments.sh --instances-dir PATH --instances-list FILE --configs-list FILE --time-limit FLOAT --num-seeds INT>0 [--percentage-random-sample INTEGER 0-100]
```

For example, to reproduce the full experiments for all MIPLIB instances and configurations tested in the article, one could run:
```
./script_reproduce_experiments.sh --instances-dir PATH_TO_MIPLIB_INSTANCE_FOLDER --instances-list "./instances/miplib.txt" --configs-list "./settings/settings_list.txt" --num-seeds 3 --time-limit 3600
```

### Random Sampling Execution

The optinal parameter [--percentage-random-sample INTEGER 0-100] determines a percentage of the instances listed in --instances-list that should be randomly selected to run. The default value is 100, case where the full instance set is considered.

### Results generation

In the last step of the script execution, the results of the execution are compiled and stored in a sub-folder within [solutions](solutions). The sub-folder is created during the execution with the current date/timestamp. Aside from a solution file for each combination of *instance x configuration x seed*, other result files are stored in this sub-folder:
- **sorted_files_list.txt**:  keeps the list of the instances considered in the experiments (specially useful when randomly sampling);
- **summary.csv**: raw results for each *instance x configuration x seed* execution. Important columns are *status* (which points out if the execution suceeded in findind a feasible solution or not) and *total time* (in seconds);
- The remaining .txt files are the *body* of the LaTeX tables displayed in the article.

## References

<a id="1">[1]</a> 
M. Fischetti, F. Glover, and A. Lodi. The feasibility pump. Mathematical
Programming, 104(1):91–104, 2005.

<a id="2">[2]</a> 
L. Bertacco, M. Fischetti, and A. Lodi. A feasibility pump heuristic for
general mixed-integer problems. Discrete Optimization, 4(1):63–76, 2007.

<a id="3">[3]</a> 
T. Achterberg and T. Berthold. Improving the feasibility pump. Discrete
Optimization,     4(1):77–86, 2007.

<a id="4">[4]</a> 
M. Fischetti and D. Salvagnin. Feasibility pump 2.0. Mathematical Programming 
Computation, 1(2-3):201–222, 2009

<a id="5">[5]</a> 
Berthold, T., Mexi, G., Salvagnin, D.: Using multiple reference vectors and objective
scaling in the feasibility pump. EURO Journal on Computational Optimization
11, 100066

<a id="6">[6]</a> 
L. Assunção, S. Urrutia, and A. C. Santos. Kernel pump: Boosting the convergence of the
feasibility pump through a kernel search-based decomposition. Mathematical Programming
Computation, 2026. doi:10.1007/s12532-026-00333-2 

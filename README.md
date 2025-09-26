# Kernel Pump and Feasibility Pump Collection

At this moment, the code supports the Feasibility Pump (FP) versions [[1]](#1),[[2]](#2),[[3]](#3),[[4]](#4),[[5]](#5) and the improved method called Kernel Pump (KP).
This code was developed by adapting and extending the freely available [FeasPumpCollection](https://github.com/GioniMexi/FeasPumpCollection).

## Code overview

All the code is available at the folder [src](src), being the file [src/main.cpp](src/main.cpp) the starting point. The main KP and FP codes are in [src/kernelpump.cpp](src/kernelpump.cpp) and [src/feaspump.cpp](src/feaspump.cpp). Interfaces to the supported LP solvers are in [src/cpxmodel.cpp](src/cpxmodel.cpp), [src/xprsmodel.cpp](src/xprsmodel.cpp), [src/scipmodel.cpp](src/scipmodel.cpp), and [src/pdlpmodel.cpp](src/pdlpmodel.cpp).
Transformers (objects responsible for rounding a solution) are in [src/transformers.cpp](src/transformers.cpp). Different strategies for sorting variables before rounding the FP are in [src/rankers.cpp](src/rankers.cpp).

## Compilation instructions

This project uses CMake. For compiling, perform the following steps on the root folder of the project. Use the correct paths on your machine for CPLEX, XPRESS, SCIP:

- mkdir build
- cd build
- cmake -DCMAKE_BUILD_TYPE=Release -S=.. -DCPLEX_ROOT_DIR=/opt/ilog/cos129/cplex -DXPRESSDIR=/opt/fico/xpressmp87 -DSCIP_DIR=/opt/scip..
- make -j12
- IF YOU NEED TO RUN IN **VERBOSE MODE**, comment line 'add_definitions(-DSILENT_EXEC)' in the main [CMAKELists.txt](CMakeLists.txt) file **BEFORE** executing 'make -j12'.

* For now, only the CPLEX interface is ready for use.

## Usage

The execution file is stored in the folder [build](build) after compiled.
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
./script_reproduce_experiments.sh --instances-dir PATH_TO_MIPLIB_INSTANCE_FOLDER --instances-list "./instances/miplip.txt" --configs-list "./settings/settings_list.txt" --num-seeds 3 --time-limit 3600
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

<a id="">[5]</a> 
Berthold, T., Mexi, G., Salvagnin, D.: Using multiple reference vectors and objective
scaling in the feasibility pump. EURO Journal on Computational Optimization
11, 100066 

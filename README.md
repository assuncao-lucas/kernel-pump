# Kernel Pump and Feasibility Pump Collection

At this moment, the code supports the Feasibility Pump versions [[1]](#1),[[2]](#2),[[3]](#3),[[4]](#4),[[5]](#5) and the improved method called Kernel Pump.
This code was developed by adapting and extending the freely available [FeasPumpCollection](https://github.com/GioniMexi/FeasPumpCollection).

Compilation instructions
------------------------

This project uses CMake. Minimal steps (use the correct paths on your machine for CPLEX, XPRESS, SCIP):

- mkdir build
- cd build
- cmake -DCMAKE_BUILD_TYPE=Release -S=.. -DCPLEX_ROOT_DIR=/opt/ilog/cos129/cplex -DXPRESSDIR=/opt/fico/xpressmp87 -DSCIP_DIR=/opt/scip..
- make -j12
- IF YOU NEED TO RUN IN SILENT MODE, uncomment line '# add_definitions(-DSILENT_EXEC)' in the main CMAKELists.txt file BEFORE executing 'make -j12'.

* For now, only the CPLEX interface is ready for use.

Code overview
-------------

The main KP and FP codes are in kernelpump.cpp and feaspump.cpp. Interfaces to the supported LP solvers are in cpxmodel.cpp, xprsmodel.cpp, scipmodel.cpp, and pdlpmodel.cpp.
Transformers (objects responsible for rounding a solution) are in transformers.cpp. Different strategies for sorting variables before rounding are in rankers.cpp.

Usage
-------------
```
$ .kp prob_file --config (-c) config_file
```

Is is also possible to overrride config arguments like this:
```
$ .kp prob_file --config (-c) config_file arg=value
```

Experiments reproduction
-------------
In the main experiments of the Kernel Pump method, three benchmark of instances were used:
- [MIPLIB 2017 benchmark](https://miplib.zib.de/downloads/benchmark.zip)
- [Robust Steiner Team Orienteering Problem instances](https://drive.google.com/file/d/1LPPfpt_mbNgHhu0PWz8pNHT4-Qvi9ACa/view?usp=drive_link)
- [Traveling Tournament Problem with Predefined Venues](https://drive.google.com/file/d/1cf_0n4XnAH7WJKVh6aOf6OI10p2Ep7c0/view?usp=drive_link)

The specific instances used from each benchmark are listed in the files:
- instances/all_problems_benchmark.txt
- instances/r-stop-dp.txt
- instances/travel-tour-problem-circ.txt

The configurations tested are stored in the 'settings/all' folder:
- "config1" corresponds to "FP<sup>*</sup>"
- "config3" corresponds to "FP<sup>-</sup>"
- "config20" corresponds to "FP$<sup>+</sup>"
- "config16" corresponds to "KP$<sup>*</sup>"
- "config17" corresponds to "KP<sup>-</sup>$"
- "config21" corresponds to "KP<sup>+</sup>$"
- "config23" corresponds to "KP<sup>+/-</sup>$"
- "config19" corresponds to "CPLEX<sub>std</sub>$"
- "config18" corresponds to "CPLEX<sub>feas</sub>$"

All the experiments were run in silent mode.
It is possible to run all the experiments by calling [this script file](script). It is necessary to change the 'instances_dir' and 'instances_list' variables of the script accordingly.

Results compilation instructions
--------------------------------

Code under folder results/

- mkdir build
- cd build
- cmake -DCMAKE_BUILD_TYPE=Release -S=.. -DCPLEX_ROOT_DIR=/opt/ilog/cos129/cplex -DXPRESSDIR=/opt/fico/xpressmp87 -DSCIP_DIR=/opt/scip..
- make -j12

References
-------------
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

R-STOP-DP/r-stop-dp-r50_4.lp * outlier
R-STOP-DP/r-stop-dp-rc100_1.lp * outlier
collection/mining.mps * outlier
collection/2club200v15p5scn.mps
collection/8div-n59k10.mps
collection/academictimetablebig.mps
collection/neos-4332801-seret.mps
collection/ivu06.mps
benchmark/radiationm40-10-02.mps
benchmark/radiationm18-12-05.mps
benchmark/neos-3656078-kumeu.mps
benchmark/s100.mps
benchmark/s250r10.mps
benchmark/physiciansched6-2.mps
benchmark/neos-5049753-cuanza.mps
benchmark/mushroom-best.mps

configs:
config1: classical FP
config2: classical FP + more randomness in rounding
config3: classical FP with objective min sum x.
config4: classical KP (try enforce LP feas initial kernel and init each bucket with less fractional solution from previous buckets)
config5: classical KP with objective min sum x. 
config6: config4 without enforcing LP feasibility at first kernel
config7: config5 without enforcing LP feasibility at first kernel

config8: config1 + 3 multiple vector.
config9: config3 + 3 multiple vector.
config10: config4 + 3 multiple vector.
config11: config5 + 3 multiple vector.
config12: config6 + 3 multiple vector.
config13: config7 + 3 multiple vector.

config14: OFP 2.0 + 3 multiple vectors.
config15: kernel pump with OFP 2.0 + 3 multiple vector.

config16: config4 with sortByFractionalPart = 1
config17: config5 with sortByFractionalPart = 1

config18: cplex with emphasis in feasibility
config19: cplex without emphasis in feasibility

config20: classical FP with objective MAX sum x.
config21: config17 with objective MAX sum x.

config22: config5 with objective MAX sum x (sortByFractionalPart = 0).

config23: config17 with objective MAX sum x for the kernel generation, BUT min sum x in the objectives of the FPs

config24: config17 with Null objective in KP and FP.

config25: config16 with reversed objective function KP and FP.
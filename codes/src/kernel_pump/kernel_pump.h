#pragma once

#include <vector>
#include <list>
#include <boost/dynamic_bitset.hpp>
#include "src/feasibility_pump/feasibility_pump.h"
#include "src/heuristic_solution.h"

class Problem;
class IloNumArray;

class KernelPump
{
public:
    explicit KernelPump() = default;
    virtual ~KernelPump();

    void Init(Problem *problem);
    bool Run(int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool sort_by_fractional_part, bool reset_fp_initial_basis_at_new_loop, bool always_force_bucket_vars_into_kernel);

private:
    FeasibilityPump feasibility_pump_;

    boost::dynamic_bitset<> curr_kernel_bitset_;
    boost::dynamic_bitset<> best_basis_;
    std::vector<boost::dynamic_bitset<>> buckets_bitsets_;

    bool found_int_feasible_solution_ = false;
    double best_basis_value_ = std::numeric_limits<double>::infinity();
    Problem *problem_ = nullptr;
    KPHeuristicSolution solution_;

    void ResetProblem();
    bool BuildKernelAndBuckets(int ks_max_size_bucket, bool sort_by_fractional_part);
    void PrintKernelAndBuckets();

    const KPHeuristicSolution &solution() const { return solution_; }
};
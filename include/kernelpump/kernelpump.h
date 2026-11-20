#pragma once

#include "kernelpump/mipmodel.h"
#include <string>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <utils/timer.h>
#include "kernelpump/feasibility_heuristic.h"

using namespace dominiqs;

class IloNumArray;

class KernelPump
{
public:
    explicit KernelPump() = default;
    virtual ~KernelPump() = default;

    void readConfig();
    bool Init(MIPModelPtr problem);
    bool Run(double time_limit);

    /** Map an fp.reOptMethod-style spelling ("default"/"primal"/"dual"/"barrier") to the solver's
     * LP method code. Throws std::runtime_error on an unknown spelling. Static so that main.cpp can
     * validate the config key up front, outside the try/catch that -DSILENT_EXEC would mute. */
    static char ParseLpMethod(const std::string &name);

    void Reset();
    bool foundSolution()
    {
        return found_int_feasible_solution_;
    }
    int getIterations() const
    {
        // main.cpp reads this even when Init() failed, i.e. before the heuristic exists;
        // the previous by-value FeasibilityPump member reported 0 in that case.
        return heuristic_ ? heuristic_->getIterations() : 0;
    }

    double getClosestDist() const
    {
        return closest_dist_;
    }

    void getClosestFrac(std::vector<double> &x) const;

    double getTimeSpentBuildingKernelBuckets() const
    {
        return time_spent_building_kernel_buckets_;
    }

    int getNumBuckets() const
    {
        return buckets_bitsets_.size();
    }

    int getLastBucketVisited() const
    {
        return last_bucket_visited_;
    }

    int getFirstBucketToIterPump() const
    {
        return first_bucket_to_iter_pump_;
    }

    int getNumVarsInKernel() const
    {
        return curr_kernel_bitset_.count();
    }

    int num_binary_vars_with_value_1_in_solution() const
    {
        return num_binary_vars_with_value_1_in_solution_;
    }

    int num_binary_vars_initial_kernel() const
    {
        return num_binary_vars_initial_kernel_;
    }

    int num_binary_vars_initial_kernel_refined() const
    {
        return num_binary_vars_initial_kernel_refined_;
    }

    void getSolution(std::vector<double> &solution) const;

private:
    // The feasibility heuristic KP drives, chosen by kp.heuristic and created at the END of a
    // successful Init(), so it exists if and only if Init() succeeded.
    // IMPORTANT: exactly one instance per KP run -- implementations may carry state across
    // buckets, and FeasibilityPump does (its iteration limits are a budget shared by all
    // buckets). See docs/spec-heuristic-interface.md section 5.
    FeasibilityHeuristicPtr heuristic_;
    boost::dynamic_bitset<> binaries_;
    boost::dynamic_bitset<> gintegers_;
    boost::dynamic_bitset<> continuous_;

    bool BuildKernelAndBuckets(double time_limit);
    int addVarToBucket(int var_index, const std::vector<double> &var_values, boost::dynamic_bitset<> &curr_bucket_bitset, boost::dynamic_bitset<> &total_added_vars_bitset) const;
    void PrintKernelAndBuckets();

    // solve data.
    MIPModelPtr model_;
    MIPModelPtr original_model_; // must be saved for converting post solve solution in case of presolve.
    StopWatch kp_watch_;
    boost::dynamic_bitset<> curr_kernel_bitset_;
    std::vector<double> closest_frac_;
    std::vector<boost::dynamic_bitset<>> buckets_bitsets_;
    int last_bucket_visited_ = 0;
    double closest_dist_ = INFBOUND;
    double time_spent_building_kernel_buckets_ = 0.0;
    double total_time_spent_ = 0.0;
    bool found_int_feasible_solution_ = INFBOUND;
    double primal_bound_ = std::numeric_limits<double>::infinity();
    bool has_presolve_ = false;
    int first_bucket_to_iter_pump_ = -1; // first bucket for which the feasibility pump is able to iterate (find problem LP feasible and not MIP infeasible in the presolve).
    std::vector<double> solution_;
    int num_binary_vars_with_value_1_in_solution_ = 0;
    int num_binary_vars_initial_kernel_ = 0;
    int num_binary_vars_initial_kernel_refined_ = 0;
    std::shared_ptr<std::vector<boost::dynamic_bitset<>>> cols_dependency_;

    // parameters.
    std::string heuristic_name_ = "feaspump"; // kp.heuristic
    // Simplex method for the LPs solved while repairing the initial kernel. Formerly taken from
    // FeasibilityPump::getReOptMethod(), which at that point in Run() had not yet had
    // readConfig() called on it and so always yielded the constructor default 'S'. Now explicit.
    // See docs/spec-heuristic-interface.md section 6.3.
    char kernel_repair_lp_method_ = 'S'; // kp.kernelRepairLpMethod
    bool try_enforce_feasibility_initial_kernel_ = false;
    bool build_kernel_based_on_null_obj_ = false;
    bool build_kernel_based_on_sum_vars_obj_ = false;
    bool build_kernel_based_on_sum_vars_obj_max_sense_ = false;
    bool reverse_obj_func_ = false;
    bool reset_fp_initial_basis_at_new_loop_ = false;
    bool sort_by_fractional_part_ = false;
    bool always_force_bucket_vars_into_kernel_ = false;
    bool buckets_by_relaxation_layers_ = false;
    bool buckets_by_variable_dependency_ = false;
    int num_bucket_layers_ = 0;
    int max_size_buckets_ = 0;
};
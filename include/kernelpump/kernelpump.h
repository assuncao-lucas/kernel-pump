#pragma once

#include "kernelpump/mipmodel.h"
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <utils/timer.h>
#include "kernelpump/feaspump.h"

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

    void Reset();
    bool foundSolution()
    {
        return found_int_feasible_solution_;
    }
    int getIterations() const
    {
        return feasibility_pump_.getIterations();
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

    void getSolution(std::vector<double> &solution) const;

private:
    FeasibilityPump feasibility_pump_;
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
    std::shared_ptr<std::vector<boost::dynamic_bitset<>>> cols_dependency_;

    // parameters.
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
#include <cmath>
#include <algorithm>
#include "kernelpump/kernelpump.h"
#include <utils/consolelog.h>
#include <utils/fileconfig.h>

using namespace dominiqs;

// DEFAULT Kernel Search parameters.
static const int K_KP_DEFAULT_NUM_BUCKET_LAYERS = 10;
static const int K_KP_DEFAULT_NAX_SIZE_BUCKETS = 100;

void KernelPump::Reset()
{
    kp_watch_.reset();
    closest_dist_ = INFBOUND;
    time_spent_building_kernel_buckets_ = 0.0;
    total_time_spent_ = 0.0;
    found_int_feasible_solution_ = false;
    primal_bound_ = INFBOUND;
    has_presolve_ = false;
    solution_.clear();
    cols_dependency_.reset();
    first_bucket_to_iter_pump_ = -1;
}

void KernelPump::getClosestFrac(std::vector<double> &x) const
{
    if (closest_frac_.empty())
    {
        x.clear();
        return;
    }
    // uncrush solution
    std::vector<double> postClosestFrac;
    if (has_presolve_)
    {
        DOMINIQS_ASSERT((int)closest_frac_.size() == model_->ncols());
        postClosestFrac = original_model_->postsolveSolution(closest_frac_);
    }
    else
        postClosestFrac = closest_frac_;

    x.resize(postClosestFrac.size());
    copy(postClosestFrac.begin(), postClosestFrac.end(), x.begin());
}

void KernelPump::readConfig()
{
    try_enforce_feasibility_initial_kernel_ = gConfig().get("kp.tryEnforceFeasibilityInitialKernel", false);
    build_kernel_based_on_null_obj_ = gConfig().get("kp.buildKernelBasedOnNullObjective", false);
    build_kernel_based_on_sum_vars_obj_ = gConfig().get("kp.buildKernelBasedOnSumVarsObjective", false);
    reverse_obj_func_ = gConfig().get("kp.reverseObjectiveFunction", false);
    build_kernel_based_on_sum_vars_obj_max_sense_ = gConfig().get("kp.buildKernelBasedOnSumVarsObjectiveMaxSense", false);
    reset_fp_initial_basis_at_new_loop_ = gConfig().get("kp.resetFPBasisAtNewPump", false);
    sort_by_fractional_part_ = gConfig().get("kp.sortByFractionalPart", false);
    always_force_bucket_vars_into_kernel_ = gConfig().get("kp.forceBucketVarsIntoKernel", false);
    buckets_by_relaxation_layers_ = gConfig().get("kp.buildBucketsByRelaxationLayers", false);
    buckets_by_variable_dependency_ = gConfig().get("kp.buildBucketsConsideringVariableDependency", false);
    num_bucket_layers_ = gConfig().get("kp.numBucketLayers", K_KP_DEFAULT_NUM_BUCKET_LAYERS);
    max_size_buckets_ = gConfig().get("kp.maxBucketSize", K_KP_DEFAULT_NAX_SIZE_BUCKETS);

    // log.
    consoleInfo("[config kp]");
    LOG_ITEM("kp.tryEnforceFeasibilityInitialKernel", try_enforce_feasibility_initial_kernel_);
    LOG_ITEM("kp.buildKernelBasedOnNullObjective", build_kernel_based_on_null_obj_);
    LOG_ITEM("kp.buildKernelBasedOnSumVarsObjective", build_kernel_based_on_sum_vars_obj_);
    LOG_ITEM("kp.reverseObjectiveFunction", reverse_obj_func_);
    LOG_ITEM("kp.buildKernelBasedOnSumVarsObjectiveMaxSense", build_kernel_based_on_sum_vars_obj_max_sense_);
    LOG_ITEM("kp.resetFPBasisAtNewPump", reset_fp_initial_basis_at_new_loop_);
    LOG_ITEM("kp.sortByFractionalPart", sort_by_fractional_part_);
    LOG_ITEM("kp.forceBucketVarsIntoKernel", always_force_bucket_vars_into_kernel_);
    LOG_ITEM("kp.buildBucketsByRelaxationLayers", buckets_by_relaxation_layers_);
    LOG_ITEM("kp.buildBucketsConsideringVariableDependency", buckets_by_variable_dependency_);
    LOG_ITEM("kp.numBucketLayers", num_bucket_layers_);
    LOG_ITEM("kp.maxBucketSize", max_size_buckets_);
}

bool KernelPump::Init(MIPModelPtr model)
{
    // INIT
    consoleInfo("[kpInit]");
    Reset();

    original_model_ = model; // must be saved for converting post solve solution in case of presolve.
    consoleLog("originalProblem: #rows={} #cols={} #nnz={}",
               original_model_->nrows(), original_model_->ncols(), original_model_->nnz());

    MIPModelPtr premodel;
    bool mipPresolve = gConfig().get("mipPresolve", true);

    if (mipPresolve)
    {
        double timeLimit = gConfig().get("timeLimit", 1e+20);
        original_model_->dblParam(DblParam::TimeLimit, timeLimit);

        if (!(original_model_->presolve())) // if fails, already halts, cause problem is infeasible!
        {
            consoleError("kpPresolvedProblem: MIP infeasible");
            return false;
        }
        premodel = original_model_->presolvedModel();
        if (!premodel)
        {
            // presolve made no reduction: just clone the original model
            consoleLog("kpPresolvedProblem: no reductions");
            premodel = original_model_->clone();
        }
        else
        {
            has_presolve_ = true;
            consoleLog("kpPresolvedProblem: #rows={} #cols={} #nnz={}",
                       premodel->nrows(), premodel->ncols(), premodel->nnz());
        }
    }
    else
    {
        // presolve disabled: just clone the original model
        premodel = original_model_->clone();
    }
    DOMINIQS_ASSERT(premodel);

    model_ = premodel;
    curr_kernel_bitset_.reset();
    buckets_bitsets_.clear();

    int num_vars = model_->ncols();
    closest_frac_.clear();
    closest_dist_ = INFBOUND;

    if (model_->objSense() == ObjSense::MIN)
        primal_bound_ = INFBOUND;
    else
        primal_bound_ = -INFBOUND;

    std::vector<char> xType(num_vars);
    model_->ctypes(&xType[0]);
    binaries_ = boost::dynamic_bitset<>(num_vars, 0);
    gintegers_ = boost::dynamic_bitset<>(num_vars, 0);
    continuous_ = boost::dynamic_bitset<>(num_vars, 0);

    // std::vector<std::string> var_names;
    // model_->colNames(var_names);

    for (int i = 0; i < num_vars; ++i)
    {
        if (xType[i] == 'B')
        {
            binaries_[i] = 1;
            // std::cout << var_names[i] << std::endl;
        }
        else if (xType[i] == 'I')
        {
            gintegers_[i] = 1;
        }
        else
            continuous_[i] = 1;
        // std::cout << var_names[i] << " " << xType[i] << std::endl;
    }

    return true;
}

void KernelPump::getSolution(std::vector<double> &x) const
{
    DOMINIQS_ASSERT(found_int_feasible_solution_);

    // uncrush solution
    std::vector<double> post_solution;
    if (has_presolve_)
    {
        DOMINIQS_ASSERT((int)solution_.size() == model_->ncols());
        post_solution = original_model_->postsolveSolution(solution_);
        original_model_->postsolve();
    }
    else
        post_solution = solution_;

    x.resize(post_solution.size());
    copy(post_solution.begin(), post_solution.end(), x.begin());
}

bool KernelPump::BuildKernelAndBuckets(double time_limit)
{
    consoleInfo("[kp build kernel/buckets]");
    if (!model_)
        return false;

    // retrieve info from original MIP model.
    int num_vars = model_->ncols();
    int num_binary_vars = binaries_.count();

    // keeps track of variables already addded to kernel or buckets.
    boost::dynamic_bitset<> total_added_vars_bitset(num_vars, 0);

    // if maximization problem, should order by NON-INCREASING order of reduced costs.
    bool invert_ordering_reduced_costs = (model_->objSense() != ObjSense::MIN);
    // if NOT sort by fractional values, should order in NON-INCREASING order of the actual relaxation values.
    bool invert_ordering_values = !sort_by_fractional_part_;

    curr_kernel_bitset_ = boost::dynamic_bitset<>(num_vars, 0);

    if (num_binary_vars == 0) // already stop building kernel if no binary var.
        return true;

    if (buckets_by_variable_dependency_)
    {
        consoleInfo("[computing vars dependency]");
        cols_dependency_ = model_->colsDependency();
        // model_->printDependencies();
    }

    // clone model and relax integrality!
    MIPModelPtr cloned_model_lp = model_->clone();

    if (build_kernel_based_on_null_obj_)
    {
        std::vector<double> coefs(num_vars, 0);
        std::vector<int> indexes(num_vars);

        for (int i = 0; i < num_vars; ++i)
            indexes[i] = i;

        cloned_model_lp->objcoefs(num_vars, &(indexes[0]), &(coefs[0]));
        cloned_model_lp->objOffset(0.0); //< get rid of offset
    }
    else if (build_kernel_based_on_sum_vars_obj_)
    {
        std::vector<double> coefs(num_vars, 0);
        std::vector<int> indexes(num_vars);

        for (int i = 0; i < num_vars; ++i)
        {
            if (binaries_[i])
                coefs[i] = 1;
            indexes[i] = i;
        }

        cloned_model_lp->objcoefs(num_vars, &(indexes[0]), &(coefs[0]));
        if (build_kernel_based_on_sum_vars_obj_max_sense_)
            cloned_model_lp->objSense(ObjSense::MAX);
        else
            cloned_model_lp->objSense(ObjSense::MIN);
        cloned_model_lp->objOffset(0.0); //< get rid of offset
    }

    // cloned_model_lp->writeModel("antes.lp");
    if (reverse_obj_func_)
    {
        const auto curr_obj_sense = cloned_model_lp->objSense();
        cloned_model_lp->objSense(curr_obj_sense == ObjSense::MAX ? ObjSense::MIN : ObjSense::MAX);
    }
    // cloned_model_lp->writeModel("depois.lp");
    // std::cout << "escreveu!" << std::endl;
    // getchar();
    // getchar();

    cloned_model_lp->switchToLP();
    cloned_model_lp->handleCtrlC(true);

    auto time_left = std::max(time_limit - kp_watch_.getElapsed(), 0.0);
    cloned_model_lp->dblParam(DblParam::TimeLimit, time_left);
    // make this initial solve with dual simplex to try to avoid finding optimal value, but with no primal solution (which might more often happen with barrier method).
    bool result = cloned_model_lp->lpopt('D', false, true);
    cloned_model_lp->handleCtrlC(false);

    bool pFeas = cloned_model_lp->isPrimalFeas();
    bool pbAborted = cloned_model_lp->aborted();
    bool pbInfeasTimeReached = cloned_model_lp->isInfeasibleOrTimeReached();

    if (pbAborted)
    {
        consoleError("kpBuild failed");
        consoleWarn("Cause: opt aborted");
        return false;
    }
    else if (!result)
    {
        consoleError("kpBuild failed");
        consoleWarn("Cause: opt failed");
        return false;
    }
    else if (pbInfeasTimeReached)
    {
        consoleError("kpBuild failed");
        consoleWarn("Cause: model infeasible or time reached");
        return false;
    }
    else if (!pFeas)
    {
        consoleError("kpBuild failed");
        consoleWarn("Cause: could not find feasible solution (but problem might be feasible)");
        return false;
    }

    // Sort (binary) vars in non-ascending order of LP values. For vars with val == 0, sort by reduced costs.
    std::vector<double> var_values(num_vars, 0);
    std::vector<double> var_reduced_costs(num_vars, 0);
    boost::dynamic_bitset<> non_zero_value_binary_vars(num_vars, 0);
    cloned_model_lp->sol(&var_values[0]);
    cloned_model_lp->reduced_costs(&var_reduced_costs[0]);

    struct VarValueReducedCost
    {
        int var_index;
        double value;
        double reduced_cost;
    };

    // only consider the binary variables in the ordering.
    std::vector<VarValueReducedCost> var_value_red_cost;
    for (int var_index = binaries_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = binaries_.find_next(var_index))
    {
        // std::cout << var_values[var_index] << std::endl;
        auto curr_value = sort_by_fractional_part_ ? fabs(round(var_values[var_index]) - var_values[var_index]) : var_values[var_index];
        var_value_red_cost.push_back(VarValueReducedCost{.var_index = var_index, .value = curr_value, .reduced_cost = var_reduced_costs[var_index]});

        // only add binary vars with non-zero values.
        if (greaterThan(curr_value, 0))
            non_zero_value_binary_vars[var_index] = 1;
    }

    // std::vector<std::string> xNames;
    // model_->colNames(xNames);
    // std::cout << "Before sorting" << std::endl;
    // for (auto &item : var_value_red_cost)
    // {
    //     std::cout << xNames[item.var_index] << " " << item.value << " " << item.reduced_cost << std::endl;
    // }

    auto compare = [invert_ordering_reduced_costs, invert_ordering_values](const VarValueReducedCost &a, const VarValueReducedCost &b)
    {
        // if ((a.vertex <= num_mandatory) && (b.vertex > num_mandatory))
        //     return true;
        // if ((a.vertex > num_mandatory) && (b.vertex <= num_mandatory))
        //     return false;
        auto coef_value = invert_ordering_values ? -1 : 1;
        auto coef_red_cost = invert_ordering_reduced_costs ? -1 : 1;

        return equal(a.value, b.value) ? coef_red_cost * a.reduced_cost < coef_red_cost * b.reduced_cost : coef_value * a.value < coef_value * b.value;
    };

    std::sort(var_value_red_cost.begin(), var_value_red_cost.end(), compare);

    // std::vector<std::string> xNames;
    // model_->colNames(xNames);
    // std::cout << "After sorting" << std::endl;
    // int num_zero_rc = 0, num_positive_rc = 0, num_negative_rc = 0;
    // for (auto &item : var_value_red_cost)
    // {
    //     if (equal(item.value, 0))
    //     {
    //         if (equal(item.reduced_cost, 0))
    //             ++num_zero_rc;
    //         if (lessThan(item.reduced_cost, 0))
    //             ++num_negative_rc;
    //         if (greaterThan(item.reduced_cost, 0))
    //             ++num_positive_rc;
    //     }
    //     std::cout << xNames[item.var_index] << " " << item.value << " " << item.reduced_cost << std::endl;
    // }
    // consoleInfo("num_zero_rc: {}, num_positive_rc: {}, num_negative_rc: {}", num_zero_rc, num_positive_rc, num_negative_rc);
    // getchar();
    // getchar();

    if (!buckets_by_relaxation_layers_)
    {
        int size_kernel = std::min(num_binary_vars, max_size_buckets_);

        for (int i = 0; i < size_kernel; ++i)
            curr_kernel_bitset_[var_value_red_cost[i].var_index] = 1;

        consoleLog("Kernel: {}/{} vars", size_kernel, num_binary_vars);

        int num_buckets = std::ceil(1.0 * ((num_binary_vars - size_kernel)) / max_size_buckets_);
        // note: even though we only consider binary variables in the problem, the bitsets keep the size of all the variables of the original problem.
        buckets_bitsets_ = std::vector<boost::dynamic_bitset<>>(num_buckets, boost::dynamic_bitset<>(num_vars, 0));

        int vars_added = size_kernel; // since already added some vertices to the kernel.
        for (int curr_bucket = 0; curr_bucket < num_buckets; ++curr_bucket)
        {
            // std::cout << "bucket " << curr_bucket << std::endl;
            int num_elements_in_bucket = std::min(max_size_buckets_, num_binary_vars - size_kernel - curr_bucket * max_size_buckets_);
            // std::cout << "num elements in bucket " << num_elements_in_bucket << std::endl;
            for (int curr_element_in_bucket = 0; curr_element_in_bucket < num_elements_in_bucket; ++curr_element_in_bucket)
            {
                buckets_bitsets_[curr_bucket][var_value_red_cost[vars_added].var_index] = 1;
                ++vars_added;
            }
            consoleLog("Bucket {}: {}/{} vars", curr_bucket, num_elements_in_bucket, num_binary_vars);
        }

        assert(vars_added == num_binary_vars);
    }
    else
    {
        buckets_bitsets_.clear();
        boost::dynamic_bitset<> curr_bucket_bitset(num_vars, 0);
        double first_value = var_value_red_cost.front().value, last_value = var_value_red_cost.back().value;
        double delta_value_per_bucket = (last_value - first_value) / num_bucket_layers_;
        int delta_sign = sign(delta_value_per_bucket); // indicates if value increases of decreases in sequence.
        consoleLog("interval: [{:.4f},{:.4f}]", first_value, last_value);
        consoleLog("delta: {:.4f}", delta_value_per_bucket);
        consoleLog("delta_sign: {}", delta_sign);

        auto start_range = first_value;
        auto end_range = first_value + delta_value_per_bucket;

        int var_count = 0;
        double curr_var_value = var_value_red_cost[var_count].value;
        int bucket_count = 0;

        // if all variables with same value, add all of them to a same bucket.
        if (equal(first_value, last_value))
        {
            int num_dependent_vars_added = 0;
            while (var_count < num_binary_vars)
            {
                int current_red_cost_sign = sign(var_value_red_cost[var_count].reduced_cost);

                while ((var_count < num_binary_vars) && (current_red_cost_sign == sign(var_value_red_cost[var_count].reduced_cost)))
                {
                    // consoleInfo("{}", var_count);
                    int num_vars_added = addVarToBucket(var_value_red_cost[var_count].var_index, var_values, curr_bucket_bitset, total_added_vars_bitset);
                    num_dependent_vars_added += std::max(num_vars_added - 1, 0); // excludes de origina var added, to count just its dependents.
                    ++var_count;
                }

                consoleLog("range: [{:.4f},{:.4f}], red cost sign: {} | num vars added: {} | {} of them added from dependency", start_range, end_range, current_red_cost_sign, curr_bucket_bitset.count(), num_dependent_vars_added);

                if (curr_bucket_bitset.count() > 0)
                {
                    if (bucket_count == 0)
                    {
                        curr_kernel_bitset_ = curr_bucket_bitset;
                        consoleLog("Kernel: {}/{} vars", curr_kernel_bitset_.count(), num_binary_vars);
                    }
                    else
                    {
                        buckets_bitsets_.push_back(curr_bucket_bitset);
                        consoleLog("Bucket {}: {}/{} vars", bucket_count, curr_bucket_bitset.count(), num_binary_vars);
                    }
                    ++bucket_count;
                }
                curr_bucket_bitset.reset();
            }
        }
        else
        {
            for (double value_layer = first_value; lessEqualThan(delta_sign * value_layer, delta_sign * last_value); value_layer += delta_value_per_bucket)
            {
                start_range = value_layer;
                end_range = start_range + delta_value_per_bucket;
                if (lessThan(delta_sign * start_range, delta_sign * last_value))
                    end_range = start_range + delta_value_per_bucket;
                else
                    end_range = delta_sign * std::numeric_limits<double>::infinity();

                int num_dependent_vars_added = 0;

                while (var_count < num_binary_vars && greaterEqualThan(delta_sign * curr_var_value, delta_sign * start_range) && lessThan(delta_sign * curr_var_value, delta_sign * end_range))
                {
                    int current_red_cost_sign = sign(var_value_red_cost[var_count].reduced_cost);
                    while (var_count < num_binary_vars && greaterEqualThan(delta_sign * curr_var_value, delta_sign * start_range) && lessThan(delta_sign * curr_var_value, delta_sign * end_range) && current_red_cost_sign == sign(var_value_red_cost[var_count].reduced_cost))
                    {
                        // consoleInfo("{}", var_count);
                        int num_vars_added = addVarToBucket(var_value_red_cost[var_count].var_index, var_values, curr_bucket_bitset, total_added_vars_bitset);
                        num_dependent_vars_added += std::max(num_vars_added - 1, 0); // excludes de original var added, to count just its dependents.
                        ++var_count;
                        if (var_count < num_binary_vars)
                            curr_var_value = var_value_red_cost[var_count].value;
                    }

                    consoleLog("range: [{:.4f},{:.4f}), red cost sign: {} | num vars added: {} | {} of them added from dependency", start_range, end_range, current_red_cost_sign, curr_bucket_bitset.count(), num_dependent_vars_added);

                    // only add bucket if not empty.
                    if (curr_bucket_bitset.count() > 0)
                    {
                        if (bucket_count == 0)
                        {
                            curr_kernel_bitset_ = curr_bucket_bitset;
                            int total_num_bin_vars_activate_for_feasibility = 0;
                            if (try_enforce_feasibility_initial_kernel_)
                            {
                                consoleInfo("[try to enforce LP feasibility to initial kernel]");
                                // // set objective to zero, to stop at first feasible solution
                                // std::vector<double> coefs(num_vars, 0);
                                // std::vector<int> indexes(num_vars);

                                // for (int i = 0; i < num_vars; ++i)
                                //     indexes[i] = i;

                                // cloned_model_lp->objcoefs(num_vars, &(indexes[0]), &(coefs[0]));

                                cloned_model_lp->handleCtrlC(true);
                                // at first, deactivate all binary variables.
                                cloned_model_lp->updateModelVarBounds(std::nullopt, binaries_);

                                boost::dynamic_bitset<> previous_kernel_bitset(num_vars, 0);
                                do
                                {
                                    // update model activating all current variables in current bucket (kernel).
                                    cloned_model_lp->updateModelVarBounds(curr_kernel_bitset_ - previous_kernel_bitset, std::nullopt);
                                    previous_kernel_bitset = curr_kernel_bitset_;

                                    auto time_left = std::max(time_limit - kp_watch_.getElapsed(), 0.0);
                                    cloned_model_lp->dblParam(DblParam::TimeLimit, time_left);
                                    result = cloned_model_lp->lpopt(feasibility_pump_.getReOptMethod(), false, true);
                                    cloned_model_lp->handleCtrlC(false);

                                    if (cloned_model_lp->aborted()) // result == false means LP is infeasible!
                                    {
                                        consoleError("kpBuild failed");
                                        return false;
                                    }

                                    // if solve successfull, check if solution found is indeed primal feasible.
                                    if (result)
                                    {
                                        result = false;
                                        if (cloned_model_lp->isPrimalFeas())
                                        {
                                            std::vector<double> sol(num_vars);
                                            cloned_model_lp->sol(&sol[0]);
                                            if (cloned_model_lp->isSolutionFeasible(sol))
                                                result = true;
                                        }
                                    }

                                    // if problem infeasible, try to fix it by activating variables that contribute to infeasibility.
                                    if (!result)
                                    {
                                        std::vector<int> conflicting_constraints;
                                        std::vector<int> conflicting_vars;
                                        time_left = std::max(time_limit - kp_watch_.getElapsed(), 0.0);

                                        if (equal(time_left, 0))
                                            break;

                                        cloned_model_lp->findSetOfConflictingVariables(non_zero_value_binary_vars - curr_kernel_bitset_, conflicting_constraints, conflicting_vars, true, time_left);

                                        time_left = std::max(time_limit - kp_watch_.getElapsed(), 0.0);

                                        if (equal(time_left, 0))
                                            break;

                                        // activate all variables in conflict
                                        int num_bin_vars_activate_for_feasibility_iter = 0;
                                        int num_dependent_vars_added_iter = 0;

                                        int num_binary_conflicting_vars = 0;
                                        for (auto conflicting_var : conflicting_vars)
                                        {
                                            if (binaries_[conflicting_var])
                                                ++num_binary_conflicting_vars;

                                            // note: a variable previously added to the kernel can be part of a conflict (even of upper bound type!)
                                            // in a later solving of the problem. So, better check if already added, just to avoid double work.
                                            if (binaries_[conflicting_var] && !(total_added_vars_bitset[conflicting_var]))
                                            {
                                                ++num_bin_vars_activate_for_feasibility_iter;
                                                int num_vars_added = addVarToBucket(conflicting_var, var_values, curr_kernel_bitset_, total_added_vars_bitset);
                                                num_dependent_vars_added_iter += std::max(num_vars_added - 1, 0); // excludes de original var added, to count just its dependents.
                                            }
                                        }
                                        // consoleWarn("{}", num_binary_conflicting_vars);

                                        num_dependent_vars_added += num_dependent_vars_added_iter;
                                        num_bin_vars_activate_for_feasibility_iter += num_dependent_vars_added_iter;
                                        total_num_bin_vars_activate_for_feasibility += num_bin_vars_activate_for_feasibility_iter;
                                        if (num_bin_vars_activate_for_feasibility_iter > 0)
                                            consoleLog(" * Added {} more vars to enforce feasibility | {} out of them added from dependency", num_bin_vars_activate_for_feasibility_iter, num_dependent_vars_added_iter);
                                    }
                                } while (!result && curr_kernel_bitset_.count() < num_binary_vars && curr_kernel_bitset_ != previous_kernel_bitset);

                                if (!result && curr_kernel_bitset_ != previous_kernel_bitset) // one last check to see if problem infeasible if added variables in the last iteration.
                                {
                                    // update model activating all current variables in current bucket (kernel).
                                    cloned_model_lp->updateModelVarBounds(curr_kernel_bitset_ - previous_kernel_bitset, std::nullopt);
                                    auto time_left = std::max(time_limit - kp_watch_.getElapsed(), 0.0);
                                    cloned_model_lp->dblParam(DblParam::TimeLimit, time_left);
                                    result = cloned_model_lp->lpopt(feasibility_pump_.getReOptMethod(), false, true);
                                }
                                // IMPORTANT NOTE: even when an LP feasible initial kernel is found, it might be the case that, when actually solving the sub-problem related to this kernel,
                                // the presolve might detect that the problem is MIP infeasible. This happens because the presolve considers the INTEGER problem. Thus, the LP feasible problem
                                // can still be MIP infeasible!
                                if (result)
                                    consoleInfo(" * Found LP feasible initial kernel");
                                else
                                    consoleWarn(" * Found LP infeasible initial kernel");

                                //     // restore type information
                                //     int n = cloned_model_lp->ncols();
                                //     std::vector<char> ctype(n, 'C');

                                //     for (int j = 0; j < n; j++)
                                //     {
                                //         if (binaries_[j])
                                //             cloned_model_lp->ctype(j, 'B');
                                //         else if (gintegers_[j])
                                //             cloned_model_lp->ctype(j, 'I');
                                //         else
                                //             cloned_model_lp->ctype(j, 'C');
                                //     }

                                //     std::vector<double> test_lbs(n);
                                //     std::vector<double> test_ubs(n);
                                //     std::vector<char> test_types(n);
                                //     cloned_model_lp->lbs(&test_lbs[0]);
                                //     cloned_model_lp->ubs(&test_ubs[0]);
                                //     cloned_model_lp->ctypes(&test_types[0]);
                                //     size_t test_num_continuous_vars = 0;
                                //     size_t test_num_integer_vars = 0;
                                //     size_t test_num_binary_vars = 0;
                                //     size_t test_num_fixed_zero_vars = 0;
                                //     size_t test_num_fixed_one_vars = 0;
                                //     size_t test_num_rows = cloned_model_lp->nrows();
                                //     for (int i = 0; i < n; i++)
                                //     {
                                //         // consoleLog("var {}",i);
                                //         // fixed variables are not considered integer variables, of any kind
                                //         // this is correct if the rounding function does not alter their value!
                                //         // this is trivially true for simple rounding, but some care must be
                                //         // taken for more elaborate strategies!!!
                                //         if (equal(test_lbs[i], test_ubs[i], 1e-6))
                                //         {
                                //             if (equal(test_lbs[i], 0, 1e-6))
                                //                 ++test_num_fixed_zero_vars;
                                //             else if (equal(test_ubs[i], 1, 1e-6))
                                //                 ++test_num_fixed_one_vars;
                                //         }
                                //         else if (test_types[i] != 'C')
                                //         {
                                //             if (test_types[i] == 'B')
                                //                 ++test_num_binary_vars;
                                //             else if (test_types[i] == 'I')
                                //                 ++test_num_integer_vars;
                                //         }
                                //         else
                                //             ++test_num_continuous_vars;
                                //     }

                                //     consoleInfo("model after build kernel: fixed 0: {} | fixed 1: {} | integer: {} | binary: {} | continuous: {} | rows: {}", test_num_fixed_zero_vars, test_num_fixed_one_vars, test_num_integer_vars, test_num_binary_vars, test_num_continuous_vars, test_num_rows);
                            }

                            consoleLog("Kernel: {}/{} vars | {} of them added to try to enforce feasibility | {} out of them added from dependency", curr_kernel_bitset_.count(), num_binary_vars, total_num_bin_vars_activate_for_feasibility, num_dependent_vars_added);
                        }
                        else
                        {
                            buckets_bitsets_.push_back(curr_bucket_bitset);
                            consoleLog("Bucket {}: {}/{} vars", bucket_count, curr_bucket_bitset.count(), num_binary_vars);
                        }

                        ++bucket_count;
                    }
                    curr_bucket_bitset.reset();

                    if (var_count >= num_binary_vars)
                        break;
                }
                if (var_count >= num_binary_vars)
                    break;
            }
        }
    }
    // getchar();
    // getchar();
    return true;
}

int KernelPump::addVarToBucket(int var_index, const std::vector<double> &var_values, boost::dynamic_bitset<> &curr_bucket_bitset, boost::dynamic_bitset<> &total_added_vars_bitset) const
{
    int num_vars_added = 0;
    if (!total_added_vars_bitset[var_index])
    {
        curr_bucket_bitset[var_index] = 1;
        total_added_vars_bitset[var_index] = 1;
        ++num_vars_added;
        // also add dependent variables if the case.
        if (buckets_by_variable_dependency_)
        {
            auto &var_dependency = (*cols_dependency_)[var_index];
            for (int dependent_var = var_dependency.find_first(); dependent_var != boost::dynamic_bitset<>::npos; dependent_var = var_dependency.find_next(dependent_var))
            {
                // only add binary variables not yet added and with relaxation value greater than zero.
                if (!total_added_vars_bitset[dependent_var] && binaries_[dependent_var] && greaterThan(var_values[dependent_var], 0))
                {
                    curr_bucket_bitset[dependent_var] = 1;
                    total_added_vars_bitset[dependent_var] = 1;
                    ++num_vars_added;
                }
            }
        }
    }
    return num_vars_added;
}

bool KernelPump::Run(double time_limit)
{
    if (!model_)
        return false;
    kp_watch_.start();

    auto num_vars = model_->ncols();
    double curr_fp_closet_dist = INFBOUND;

    // build Kernel by solving LP of given problem.
    double time_left = std::max(time_limit - kp_watch_.getElapsed(), 0.0);
    bool builtKernel = BuildKernelAndBuckets(time_left);
    time_spent_building_kernel_buckets_ = kp_watch_.getElapsed();
    if (builtKernel)
    {
        // getchar();
        // getchar();
        // PrintKernelAndBuckets();

        // init feasibility pump.
        feasibility_pump_.readConfig();

        // initially, deactivate all the binary variables.
        model_->updateModelVarBounds(std::nullopt, binaries_);

        int curr_bucket_index = -1; // starts from kernel.
        int total_num_buckets = buckets_bitsets_.size();
        time_left = std::max(time_limit - kp_watch_.getElapsed(), 0.0);
        double min_time_per_bucket = time_left / (total_num_buckets + 1); // +1 to consider first kernel too.
        double curr_time_limit_iteration = 0.0;

        auto curr_reference_kernel = curr_kernel_bitset_;
        auto curr_vars_entering_kernel = curr_kernel_bitset_;
        auto curr_vars_leaving_reference_kernel = boost::dynamic_bitset<>(num_vars, 0);

        for (curr_bucket_index = -1; curr_bucket_index < total_num_buckets; ++curr_bucket_index)
        {
            // if last bucket, leave running for all the remaining time left.
            curr_time_limit_iteration = (curr_bucket_index == total_num_buckets - 1) ? std::max(time_limit - kp_watch_.getElapsed(), 0.0) : min_time_per_bucket;

            if (model_->aborted() || lessEqualThan(curr_time_limit_iteration, 0.0))
                break;

            // std::cout << "curr_time_limit_iteration: " << curr_time_limit_iteration << std::endl;
            // update the reference kernel to the current kernel (+ current bucket, if not the first iteration).
            if (curr_bucket_index >= 0)
            {
                curr_reference_kernel = curr_kernel_bitset_ | buckets_bitsets_[curr_bucket_index];
                curr_vars_entering_kernel = buckets_bitsets_[curr_bucket_index];
            }

            // std::cout << " current bucket ";
            // curr_bucket_index >= 0 ? std::cout << buckets_bitsets_[curr_bucket_index] << std::endl : std::cout << " - " << std::endl;
            // std::cout << " kernel: " << curr_kernel_bitset_ << std::endl;
            // std::cout << " reference kernel: " << curr_reference_kernel << std::endl;
            // std::cout << " IN kernel: " << curr_vars_entering_kernel << std::endl;
            // std::cout << " OUT ref kernel: " << curr_vars_leaving_reference_kernel << std::endl;
            // std::cout << " curr sol: " << curr_int_y_ << std::endl;
            // std::cout << " best cost: " << curr_best_solution_value_ << std::endl;
            // cplex_.exportModel("model_before.lp");
            //  Enable in the model the variables that are active in the Kernel.

            model_->updateModelVarBounds(curr_vars_entering_kernel, curr_vars_leaving_reference_kernel);
            // std::cout << " bucket index " << curr_bucket_index << " num active bin vars: " << model_->curr_active_binary_vars().count() << std::endl;

            // std::cout << "current reference kernel: " << std::endl;
            // for (int i = 0; i < num_vars; ++i)
            // {
            //     if (curr_reference_kernel[i])
            //         std::cout << (model_->vars())[i] << " ";
            // }
            // std::cout << std::endl;
            // getchar();
            // getchar();
            // Run current Feasibility pump subproblem.
            // for each feasibility pump subproblem, retrieve the best new basis bound.
            // Note: if all iterations of the subproblem lead to infeasible basis, set bound as infinity. If at least one is solved to optimality, retrieve
            // the best bound among them.
            // Note: should never be unbounded in case of standard feasibility pump, because is simply a minimization of distance function (so, at least zero)
            // BUT, if objective feasibility pump, it might be undounded due to the influence of the original objective and it's non binary variables.

            // should continue even if feasibility pump fails (or is infeasible).
            // recall that an infeasible subproblem does not imply in infeasible original problem.
            // BUT should add all the current bucket's variables to the kernel, as to avoid later infeasibility.

            if (curr_bucket_index == -1)
                consoleInfo("[kp initial kernel]");
            else
                consoleInfo("[Kp bucket {}/{}]", curr_bucket_index + 1, total_num_buckets);

            consoleLog("#active bin vars : {}/{}", curr_reference_kernel.count(), binaries_.count());
            bool found_int_feasible_solution = false;
            bool feasible_fp = false;
            // if init fails, means that the proble is already infeasible for the current bucket.
            if (feasibility_pump_.init(model_)) // at fp.init, the original MIP problem is possibly presolved and linearly relaxed into an LP.
            {
                // for the last bucket, disable stopping criteria of FP related to max number of iterations without at least x% improvements.
                bool stop_with_no_impr_limit = (curr_bucket_index != total_num_buckets - 1);
                // give warm start to FP (with best frac basis found so far) only if a basis exists AND reset_fp_initial_basis_at_new_loop_ == false.
                std::vector<double> xStartPoint, xStartFrac;
                double xStartDist = INFBOUND;
                bool lp_primal_feas = false;
                if (!reset_fp_initial_basis_at_new_loop_ && !closest_frac_.empty())
                {
                    xStartFrac = closest_frac_;
                    xStartDist = closest_dist_;
                    lp_primal_feas = true;
                }

                StopWatch tp1;
                tp1.start();

                auto result = feasibility_pump_.pump(curr_time_limit_iteration, stop_with_no_impr_limit, xStartFrac, xStartDist, lp_primal_feas);

                // if (greaterThan(tp1.getElapsed(), curr_time_limit_iteration))
                // std::cout << " * KP  bckt " << curr_bucket_index << " " << tp1.getElapsed() << " x " << curr_time_limit_iteration << std::endl;

                found_int_feasible_solution = std::get<0>(result);
                feasible_fp = std::get<1>(result);
            }

            if (feasible_fp && first_bucket_to_iter_pump_ == -1)
                first_bucket_to_iter_pump_ = curr_bucket_index + 1;

            if (found_int_feasible_solution)
            {
                found_int_feasible_solution_ = true;
                feasibility_pump_.getSolution(solution_);
                primal_bound_ = feasibility_pump_.getPrimalBound();
                closest_dist_ = feasibility_pump_.getClosestDist();
                ++curr_bucket_index;
                curr_kernel_bitset_ = curr_reference_kernel;
                num_binary_vars_with_value_1_in_solution_ = 0;
                for (int var_index = binaries_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = binaries_.find_next(var_index))
                {
                    if (equal(solution_[var_index], 1.0))
                        ++num_binary_vars_with_value_1_in_solution_;
                }
                // std::cout << "!!!!!!!!!!!!!!!!! " << curr_kernel_bitset_.count() << std::endl;
                break;
            }
            else
            {
                if (!feasible_fp)
                {
                    // curr_vars_entering_kernel = curr_reference_kernel - curr_kernel_bitset_;
                    curr_kernel_bitset_ = curr_reference_kernel;
                    curr_vars_leaving_reference_kernel.reset();
                }
                else
                {
                    bool found_new_closest_point = false;
                    // if found a new better basis (i.e., a basis with smaller dist from integer), update best basis.
                    curr_fp_closet_dist = feasibility_pump_.getClosestDist();
                    // std::cout << curr_fp_closet_dist << " x " << closest_dist_ << std::endl;

                    if (lessThan(curr_fp_closet_dist, closest_dist_))
                    {

                        closest_dist_ = curr_fp_closet_dist;
                        feasibility_pump_.getClosestFrac(closest_frac_);
                        found_new_closest_point = true;
                    }

                    // update kernel!
                    // Note: if always_force_bucket_vars_into_kernel == true, also add all the bucket's variables to kernel.
                    if (always_force_bucket_vars_into_kernel_)
                    {
                        // curr_vars_entering_kernel = curr_reference_kernel - curr_kernel_bitset_;
                        curr_kernel_bitset_ = curr_reference_kernel;
                        curr_vars_leaving_reference_kernel.reset();
                    }
                    else if (found_new_closest_point)
                    {
                        // std::cout << "current best basis: " << std::endl;
                        // for (int i = 0; i < num_vars; ++i)
                        // {
                        //     if (closest_point_[i])
                        //         std::cout << (model_->vars())[i] << " ";
                        // }
                        // std::cout << std::endl;
                        // curr_vars_entering_kernel = closest_point_ - curr_kernel_bitset_;

                        boost::dynamic_bitset<> closest_point_bitset(num_vars, 0);
                        for (int var_index = binaries_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = binaries_.find_next(var_index))
                        {
                            if (greaterThan(closest_frac_[var_index], 0))
                                closest_point_bitset[var_index] = 1;
                        }

                        curr_kernel_bitset_ |= closest_point_bitset;

                        curr_vars_leaving_reference_kernel = curr_reference_kernel - curr_kernel_bitset_;
                    }
                    else // if hasn't found a batter basis, nothing changes in the kernel...we only remove from reference kernel the vars added in this iteration.
                    {
                        // curr_vars_entering_kernel.reset();

                        curr_vars_leaving_reference_kernel = curr_reference_kernel - curr_kernel_bitset_;
                    }
                }
            }
        }

        last_bucket_visited_ = curr_bucket_index;
        consoleLog("");
        consoleInfo("[kp results]");
        consoleLog("primalBound = {}", primal_bound_);
        consoleLog("numSols = {}", (int)found_int_feasible_solution_);
        consoleLog("lastBucketVisited = {}/{} (original kernel index == 0)", curr_bucket_index, total_num_buckets);
        consoleLog("firstBucketToIterPump = {}", first_bucket_to_iter_pump_);
        consoleLog("buildKernelAndBucketsTime = {}", time_spent_building_kernel_buckets_);
        consoleLog("totalTime = {}", kp_watch_.getElapsed());
        // std::cout << "terminou com tempo " << kp_watch_.getElapsed();
    }

    // std::cout << "Best basis found: " << curr_closest_dist_ << " " << closest_point_ << std::endl;
    // std::cout << "Elapsed time: " << timer->CurrentElapsedTime(ti) << std::endl;

    kp_watch_.stop();
    total_time_spent_ = kp_watch_.getTotal();

    // std::cout << *solution << std::endl;

    return true;
}

void KernelPump::PrintKernelAndBuckets()
{
    std::vector<std::string> xNames;
    model_->colNames(xNames);

    std::cout << "Kernel: ";

    for (int var_index = curr_kernel_bitset_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = curr_kernel_bitset_.find_next(var_index))
        std::cout << xNames[var_index] << " ";
    std::cout << std::endl;

    for (int j = 0; j < buckets_bitsets_.size(); ++j)
    {
        std::cout << "Bucket " << j + 1 << ": ";
        for (int var_index = buckets_bitsets_[j].find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = buckets_bitsets_[j].find_next(var_index))
            std::cout << xNames[var_index] << " ";

        std::cout << std::endl;
    }
}
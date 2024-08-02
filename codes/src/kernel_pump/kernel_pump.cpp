#include <cmath>
#include <algorithm>
#include "src/kernel_pump/kernel_pump.h"
#include "src/general.h"
#include "src/timer.h"
#include "src/problem.h"

KernelPump::~KernelPump()
{
    ResetProblem();
}

void KernelPump::Init(Problem *problem)
{
    Timestamp *ti = NewTimestamp();
    Timer *timer = GetTimer();
    timer->Clock(ti);

    problem_ = problem;
    auto env = problem_->env();

    curr_kernel_bitset_.reset();
    buckets_bitsets_.clear();

    bool found_int_feasible_solution_ = false;
    double best_basis_value_ = std::numeric_limits<double>::infinity();
    auto num_vars = problem_->num_vars();
    best_basis_ = boost::dynamic_bitset<>(num_vars, 0);
    best_basis_value_ = std::numeric_limits<double>::infinity();

    solution_.Reset(num_vars);

    // init feasibility pump
    feasibility_pump_.Init(problem_);

    delete ti;
    ti = nullptr;
}

void KernelPump::ResetProblem()
{
    // no need to delete, because this is an outside independent object.
    problem_ = nullptr;
}

bool KernelPump::BuildKernelAndBuckets(int ks_max_size_bucket, bool sort_by_fractional_part)
{
    if (problem_ == nullptr)
        return false;

    auto num_vars = problem_->num_vars();
    auto num_binary_vars = problem_->num_binary_vars();

    curr_kernel_bitset_ = boost::dynamic_bitset<>(num_vars, 0);

    // if maximization problem, should order by NON-INCREASING order of reduced costs.
    bool invert_ordering_reduced_costs = !(problem_->IsMinimization());
    // if NOT sort by fractional values, should order in NON-INCREASING order of the actual relaxation values.
    bool invert_ordering_values = !sort_by_fractional_part;

    problem_->Solve(true, false);

    // Sort vertices in non-ascending order of y values. For vertices with y == 0, sort by reduced costs.
    auto env = problem_->env();
    IloNumArray var_values(*env);
    IloNumArray var_reduced_costs(*env);

    problem_->GetValues(var_values);
    problem_->GetReducedCosts(var_reduced_costs);

    struct VarValueReducedCost
    {
        int var_index;
        double value;
        double reduced_cost;
    };

    auto binary_vars = problem_->binary_vars();
    // only consider the binary variables in the ordering.
    std::vector<VarValueReducedCost> vertex_value_red_cost;
    for (int var_index = binary_vars.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = binary_vars.find_next(var_index))
    {
        // std::cout << var_values[var_index] << std::endl;
        auto curr_value = sort_by_fractional_part ? fabs(round(var_values[var_index]) - var_values[var_index]) : var_values[var_index];
        vertex_value_red_cost.push_back(VarValueReducedCost{.var_index = var_index, .value = curr_value, .reduced_cost = var_reduced_costs[var_index]});
    }

    // free memory.
    var_values.end();
    var_reduced_costs.end();

    std::cout << "Before sorting" << std::endl;
    for (auto &item : vertex_value_red_cost)
    {
        std::cout << item.var_index << " " << item.value << " " << item.reduced_cost << std::endl;
    }

    auto compare = [invert_ordering_reduced_costs, invert_ordering_values](const VarValueReducedCost &a, const VarValueReducedCost &b)
    {
        // if ((a.vertex <= num_mandatory) && (b.vertex > num_mandatory))
        //     return true;
        // if ((a.vertex > num_mandatory) && (b.vertex <= num_mandatory))
        //     return false;
        auto coef_value = invert_ordering_values ? -1 : 1;
        auto coef_red_cost = invert_ordering_reduced_costs ? -1 : 1;

        return double_equals(a.value, b.value) ? coef_red_cost * a.reduced_cost < coef_red_cost * b.reduced_cost : coef_value * a.value < coef_value * b.value;
    };

    std::sort(vertex_value_red_cost.begin(), vertex_value_red_cost.end(), compare);

    std::cout << "After sorting" << std::endl;
    for (auto &item : vertex_value_red_cost)
    {
        std::cout << item.var_index << " " << item.value << " " << item.reduced_cost << std::endl;
    }

    // force that all the variables with best possible value (regardless of the ordering criteria) to be part of the initial kernel.
    int num_vars_best_value = 0;
    double best_value = -1;
    for (auto &item : vertex_value_red_cost)
    {
        if (num_vars_best_value == 0)
            best_value = item.value;

        if (double_equals(best_value, item.value))
            ++num_vars_best_value;
        else
            break;
    }

    std::cout << " **** num_vars_best_value: " << num_vars_best_value << std::endl;
    // getchar();
    // getchar();
    // start by building kernel with at least all the mandatory vertices (+ origin).
    int size_kernel = std::min(num_binary_vars, std::max(ks_max_size_bucket, num_vars_best_value));

    for (int i = 0; i < size_kernel; ++i)
        curr_kernel_bitset_[vertex_value_red_cost[i].var_index] = 1;

    int num_buckets = std::ceil(1.0 * ((num_binary_vars - size_kernel)) / ks_max_size_bucket);
    // note: even though we only consider binary variables in the problem, the bitsets keep the size of all the variables of the original problem.
    buckets_bitsets_ = std::vector<boost::dynamic_bitset<>>(num_buckets, boost::dynamic_bitset<>(num_vars, 0));

    int vars_added = size_kernel; // since already added some vertices to the kernel.
    for (int curr_bucket = 0; curr_bucket < num_buckets; ++curr_bucket)
    {
        // std::cout << "bucket " << curr_bucket << std::endl;
        int num_elements_in_bucket = std::min(ks_max_size_bucket, num_binary_vars - size_kernel - curr_bucket * ks_max_size_bucket);
        // std::cout << "num elements in bucket " << num_elements_in_bucket << std::endl;
        for (int curr_element_in_bucket = 0; curr_element_in_bucket < num_elements_in_bucket; ++curr_element_in_bucket)
        {
            buckets_bitsets_[curr_bucket][vertex_value_red_cost[vars_added].var_index] = 1;
            ++vars_added;
        }
    }

    assert(vars_added == num_binary_vars);
    return true;
}

bool KernelPump::Run(int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor, bool sort_by_fractional_part, bool reset_fp_initial_basis_at_new_loop, bool always_force_bucket_vars_into_kernel)
{
    if (problem_ == nullptr)
        return false;
    Timestamp *ti = NewTimestamp();
    Timer *timer = GetTimer();
    timer->Clock(ti);

    auto num_vars = problem_->num_vars();
    double curr_fp_sol_value = std::numeric_limits<double>::infinity();

    // build Kernel by solving LP of given problem.

    if (BuildKernelAndBuckets(ks_max_size_bucket, sort_by_fractional_part))
    {
        solution_.time_spent_building_kernel_buckets_ += timer->CurrentElapsedTime(ti);
        PrintKernelAndBuckets();

        // initially, deactivate all the binary variables.
        problem_->DeactivateAllBinaryVariables();

        double curr_time_limit_iteration = ks_max_time_limit;

        int curr_bucket_index = -1; // starts from kernel.
        int total_num_buckets = buckets_bitsets_.size();

        auto curr_reference_kernel = curr_kernel_bitset_;
        auto curr_vars_entering_kernel = curr_kernel_bitset_;
        auto curr_vars_leaving_reference_kernel = boost::dynamic_bitset<>(num_vars, 0);

        for (int curr_bucket_index = -1; curr_bucket_index < total_num_buckets; ++curr_bucket_index)
        {
            // std::cout << "curr_time_limit_iteration: " << curr_time_limit_iteration << std::endl;
            // update the reference kernel to the current kernel (+ current bucket, if not the first iteration).
            if (curr_bucket_index >= 0)
            {
                curr_reference_kernel = curr_kernel_bitset_ | buckets_bitsets_[curr_bucket_index];
                curr_vars_entering_kernel |= buckets_bitsets_[curr_bucket_index];
            }

            std::cout << " bucket index " << curr_bucket_index << std::endl;
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

            problem_->UpdateModelVarBounds(curr_vars_entering_kernel, curr_vars_leaving_reference_kernel);

            // Run current Feasibility pump subproblem.
            // for each feasibility pump subproblem, retrieve the best new basis bound.
            // Note: if all iterations of the subproblem lead to infeasible basis, set bound as infinity. If at least one is solved to optimality, retrieve
            // the best bound among them.
            // Note: should never be unbounded in case of standard feasibility pump, because is simply a minimization of distance function (so, at least zero)
            // BUT, if objective feasibility pump, it might be undounded due to the influence of the original objective and it's non binary variables.

            // should continue even if feasibility pump fails (or is infeasible).
            // recall that an infeasible subproblem does not imply in infeasible original problem.
            // BUT should add all the current bucket's variables to the kernel, as to avoid later infeasibility.
            auto feasible_fp = feasibility_pump_.Run(reset_fp_initial_basis_at_new_loop, curr_time_limit_iteration);
            if (!feasible_fp)
            {
                curr_vars_entering_kernel = curr_reference_kernel - curr_kernel_bitset_;
                curr_kernel_bitset_ = curr_reference_kernel;
                curr_vars_leaving_reference_kernel.reset();
            }
            else
            {
                // Note: if always_force_bucket_vars_into_kernel == true, also add all the bucket's variables to kernel.
                if (always_force_bucket_vars_into_kernel)
                {
                    curr_vars_entering_kernel = curr_reference_kernel - curr_kernel_bitset_;
                    curr_kernel_bitset_ = curr_reference_kernel;
                    curr_vars_leaving_reference_kernel.reset();
                }
                // if found a new better basis (i.e., a basis with smaller normalized integrality gap),
                // update kernel with the possibly new variables used in the new basis found.
                curr_fp_sol_value = feasibility_pump_.best_normalized_integrality_gap();
                std::cout << curr_fp_sol_value << " x " << best_basis_value_ << std::endl;
                if (double_less(curr_fp_sol_value, best_basis_value_))
                {
                    best_basis_value_ = curr_fp_sol_value;
                    best_basis_ = feasibility_pump_.best_basis();

                    curr_vars_entering_kernel = best_basis_ - curr_kernel_bitset_;
                    curr_kernel_bitset_ |= curr_vars_entering_kernel;

                    curr_vars_leaving_reference_kernel = curr_reference_kernel - curr_kernel_bitset_;
                }
                else // if hasn't found a batter basis, nothing changes in the kernel...we only remove from reference kernel the vars added in this iteration.
                {
                    curr_vars_entering_kernel.reset();

                    curr_vars_leaving_reference_kernel = curr_reference_kernel - curr_kernel_bitset_;
                }
            }

            if (feasibility_pump_.found_int_basis())
            {
                std::cout << "ACHOU, MISERAVEL! bucket " << curr_bucket_index << std::endl;
                solution_.is_feasible_ = true;
                solution_.found_integer_ = true;
                break;
                // BuildHeuristicSolution(solution);
            }

            curr_time_limit_iteration = std::max(curr_time_limit_iteration * ks_decay_factor, 1.0 * ks_min_time_limit);
        }
    }

    // std::cout << "Best basis found: " << curr_best_basis_value_ << " " << best_basis_ << std::endl;
    // std::cout << "Elapsed time: " << timer->CurrentElapsedTime(ti) << std::endl;

    solution_.total_time_spent_ = timer->CurrentElapsedTime(ti);

    // std::cout << *solution << std::endl;

    delete (ti);
    ti = nullptr;

    return true;
}

void KernelPump::PrintKernelAndBuckets()
{
    std::cout << "Kernel: ";

    for (int var_index = curr_kernel_bitset_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = curr_kernel_bitset_.find_next(var_index))
        std::cout << var_index << " ";
    std::cout << std::endl;

    for (int j = 0; j < buckets_bitsets_.size(); ++j)
    {
        std::cout << "Bucket " << j << ": ";
        for (int var_index = buckets_bitsets_[j].find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = buckets_bitsets_[j].find_next(var_index))
            std::cout << var_index << " ";

        std::cout << std::endl;
    }
}
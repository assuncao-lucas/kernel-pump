#include <cmath>
#include <algorithm>
#include "src/feasibility_pump/feasibility_pump.h"
#include "src/general.h"
#include "src/timer.h"

bool compare_func3(const std::pair<int, double> &v1, const std::pair<int, double> &v2)
{
	return double_greater(v1.second, v2.second);
}

FeasibilityPump::~FeasibilityPump()
{
	Reset();
}

void FeasibilityPump::Reset()
{
	if (curr_relax_basis_)
	{
		// ATTENTION! should leave it to the Problem object dealocate!
		// curr_relax_basis_->end();
		delete curr_relax_basis_;
		curr_relax_basis_ = nullptr;
	}
}

void FeasibilityPump::Init(Problem *problem)
{
	Reset();
	problem_ = problem;

	const auto num_vars = problem_->num_vars();
	solution_.Reset(num_vars);

	curr_alpha_ = 0.0;

	if (curr_relax_basis_)
	{
		curr_relax_basis_->end();
		delete curr_relax_basis_;
		curr_relax_basis_ = nullptr;
	}

	curr_relax_basis_ = new IloNumArray(*(problem_->env()));
	curr_int_basis_ = boost::dynamic_bitset<>(num_vars, 0);
	best_basis_ = boost::dynamic_bitset<>(num_vars, 0);
}

void FeasibilityPump::RetrieveAndRoundBinaryVarsValues()
{
	problem_->GetValues(*curr_relax_basis_);

	found_int_basis_ = true;
	previous_normalized_integrality_gap_ = curr_normalized_integrality_gap_;
	curr_normalized_integrality_gap_ = 0.0;

	previous_int_basis_ = curr_int_basis_;
	(curr_int_basis_).reset();
	(curr_integrality_gaps_).clear();
	auto curr_complete_int_basis = curr_int_basis_;

	// only need to check integrality of the originally binary variables.
	const auto &active_binary_vars = problem_->curr_active_binary_vars();
	for (int var_index = active_binary_vars.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = active_binary_vars.find_next(var_index))
	{
		// Note: if (*curr_relax_basis_)[var_index] == 0, there is no need to take into account this variable in the computation of the integrality gap,
		// even if the reference integer curr_int_basis_[var_index] == 1. Because there is no actual gap of integrality for this variable (it is zero),
		// even if the reference integer would try to bring it to 1 instead.
		if (!double_equals((*curr_relax_basis_)[var_index], 0.0))
		{
			curr_complete_int_basis[var_index] = 1;
			double curr_integrality_gap_var = 0.0;
			if (round((*curr_relax_basis_)[var_index])) // round == 1.
				(curr_int_basis_)[var_index] = 1;

			if (!double_equals((*curr_relax_basis_)[var_index], (curr_int_basis_)[var_index]))
			{
				found_int_basis_ = false;
				// std::cout << "var " << var_index << " relax: " << (*curr_relax_basis_)[var_index] << " round: " << (curr_int_basis_)[var_index] << std::endl;
				curr_integrality_gap_var = fabs((*curr_relax_basis_)[var_index] - 1.0 * ((curr_int_basis_)[var_index]));
				curr_normalized_integrality_gap_ += curr_integrality_gap_var;
			}
			curr_integrality_gaps_.push_back(std::pair<int, double>(var_index, curr_integrality_gap_var));
		}
	}

	if (found_int_basis_)
	{
		std::cout << "INTEIRA!!" << std::endl;
		// for (int var_index = active_binary_vars.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = active_binary_vars.find_next(var_index))
		// {
		// 	std::cout << var_index << " " << (*curr_relax_basis_)[var_index] << std::endl;
		// }
	}
	// getchar();
	// getchar();
	// normalize the expression x_1 + x_2 + ... + x_n, where n == number of ACTIVE binary variables.
	// normalization is necessary to be more fair when comparing the basis obtained for different number os active variables.
	curr_normalized_integrality_gap_ /= sqrt(problem_->curr_active_binary_vars().count());

	// update best basis if found one.
	if (double_less(curr_normalized_integrality_gap_, best_normalized_integrality_gap_))
	{
		best_normalized_integrality_gap_ = curr_normalized_integrality_gap_;
		best_basis_ = curr_complete_int_basis;
	}

	// save the basis with best known normalized integrality gap so far.
}

void FeasibilityPump::SetNewObjStage()
{
	IloExpr new_obj(*(problem_->env()));
	const auto active_binary_vars = problem_->curr_active_binary_vars();
	auto &vars = problem_->vars();

	for (int var_index = active_binary_vars.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = active_binary_vars.find_next(var_index))
		(curr_int_basis_)[var_index] == 0 ? new_obj += vars[var_index] : new_obj += (1 - vars[var_index]);

	// std::cout << new_obj << std::endl;
	// getchar();
	// getchar();
	// problem_->SetObjectiveExpression(operator*((1.0 - curr_alpha_), new_obj));
	problem_->SetObjectiveExpression(operator*((1.0 - curr_alpha_) / sqrt(problem_->num_binary_vars()), new_obj) - operator*((curr_alpha_) / problem_->original_obj_norm(), problem_->original_obj_expr()));
	new_obj.end();
}

bool FeasibilityPump::Run(bool reset_fp_initial_basis_at_new_loop, double time_limit)
{
	if (problem_ == nullptr)
		return false;
	Timestamp *ti = NewTimestamp();
	Timer *timer = GetTimer();
	timer->Clock(ti);

	best_normalized_integrality_gap_ = std::numeric_limits<double>::infinity();
	previous_normalized_integrality_gap_ = std::numeric_limits<double>::infinity();
	curr_normalized_integrality_gap_ = std::numeric_limits<double>::infinity();

	int stage_0_iter = 1, stage_2_iter = 0;
	int num_binary_vars = problem_->num_binary_vars();

	int num_flips_basis = 0, curr_num_flips = 0;
	if (K_FLIP_BASIS)
		num_flips_basis = K_FLIP_BASIS;
	else
		num_flips_basis = ceil(perturbation_flip_percentage * num_binary_vars);
	int curr_prob_of_flipping = 0;
	int num_perturbations_stage1 = 0, num_perturbations_stage2 = 0, num_restarts_stage1 = 0, num_restarts_stage2 = 0;
	(curr_alpha_) = initial_alpha_stage2;

	curr_int_basis_ = best_basis_;
	if (reset_fp_initial_basis_at_new_loop)
		curr_int_basis_.reset();

	// Stage 0: solve LP relaxation with original objective.
	if (curr_int_basis_.count() == 0) // if no initial integer basis is given, run LP with original objective function.
	{
		problem_->SetObjectiveExpression(problem_->original_obj_expr());
		problem_->IsMinimization() ? problem_->SetObjectiveSense(IloObjective::Sense::Minimize) : problem_->SetObjectiveSense(IloObjective::Sense::Maximize);
	}
	else // set minimization objective based on curr_int_basis_.
	{
		problem_->SetObjectiveSense(IloObjective::Sense::Minimize);
		SetNewObjStage();
	}

	// check if it is running the original problem (i.e., with all binary variables active).
	bool is_running_original_problem_with_all_variables_active = problem_->curr_active_binary_vars().count() == num_binary_vars;

	// std::cout << "before solve" << std::endl;
	if (!(problem_->Solve(true, false)))
	{
		// We cannot state that the original problem is infeasible with we are running just a subproblem (with some variables inactive).
		if (is_running_original_problem_with_all_variables_active)
		{
			auto curr_cplex_status = problem_->curr_status();
			if (curr_cplex_status == IloCplex::Infeasible || curr_cplex_status == IloCplex::InfOrUnbd)
				(solution_).is_infeasible_ = true;
		}
		(solution_).time_stage2_ = timer->CurrentElapsedTime(ti);
		delete ti;
		ti = nullptr;
		std::cout << "failed FP" << std::endl;
		std::cout << "num active bin vars: " << problem_->curr_active_binary_vars().count() << std::endl;
		return false;
	}
	// std::cout << "after successful solve" << std::endl;

	RetrieveAndRoundBinaryVarsValues();

	if (!(found_int_basis_) && !double_greater(timer->CurrentElapsedTime(ti), time_limit))
	{
		// Make sure that the objective sense is from now on set to minimization.
		problem_->SetObjectiveSense(IloObjective::Minimize);
		do
		{
			previous_alpha_ = curr_alpha_;
			(curr_alpha_) *= alpha_decrease_rate;
			// Update objective function according to distance from current rounded (integer) values.
			SetNewObjStage();

			++stage_2_iter;
			// std::cout << "iter: " << stage_2_iter << std::endl;

			// std::cout << "alpha: " << curr_alpha_ << std::endl;

			//  resolve updated model.
			// std::cout << "before resolve" << std::endl;
			if (!(problem_->Solve(true, false)))
			{
				// We cannot state that the original problem is infeasible with we are running just a subproblem (with some variables inactive).
				if (is_running_original_problem_with_all_variables_active)
				{
					auto curr_cplex_status = problem_->curr_status();
					if (curr_cplex_status == IloCplex::Infeasible || curr_cplex_status == IloCplex::InfOrUnbd)
						(solution_).is_infeasible_ = true;

					(solution_).time_stage2_ = timer->CurrentElapsedTime(ti);

					(solution_).num_iterations_stage2_ = stage_2_iter;
					(solution_).num_perturbations_stage2_ = num_perturbations_stage2;
					(solution_).num_restarts_stage2_ = num_restarts_stage2;
					delete ti;
					ti = nullptr;
					std::cout << "failed FP" << std::endl;
					std::cout << "num active bin vars: " << problem_->curr_active_binary_vars().count() << std::endl;
					return false;
				}
			}
			// std::cout << "after successful resolve" << std::endl;

			// update rounded values.
			RetrieveAndRoundBinaryVarsValues();

			if (!(found_int_basis_))
			{
				if ((stage_2_iter > 1) && (double_less(fabs(curr_alpha_ - previous_alpha_), K_APLHA_DECREMENT_PRECISION)) && (curr_int_basis_ == previous_int_basis_))
				{
					// std::cout << "perturbation!" << std::endl;
					++num_perturbations_stage2;

					// for (auto it = (curr_integrality_gaps_).begin(); it != (curr_integrality_gaps_).end(); ++it)
					// {
					// 	auto pos = it->first;
					// 	// curr_prob_of_flipping = it->second * 100.0; // the probability of flipping is proportional to how far the value is from an integer bound.
					// 	// std::cout << "flip prob: " << curr_prob_of_flipping << std::endl;
					// 	// int curr_rand = rand()%101;
					// 	// std::cout << "curr rand: " << curr_rand << std::endl;
					// 	// if (rand() % 101 <= curr_prob_of_flipping)
					const auto &active_binary_vars = problem_->curr_active_binary_vars();
					for (int var_index = active_binary_vars.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = active_binary_vars.find_next(var_index))
					{
						if (rand() % 2 == 1)
						{
							// std::cout << "flippou " << arc_pos << std::endl;
							((curr_int_basis_)[var_index]).flip();
						}
					}
				}
			}
		} while ((!(found_int_basis_)) && (stage_2_iter < max_iter_stage2) && !double_greater(timer->CurrentElapsedTime(ti), time_limit));
	}

	if (found_int_basis_)
	{
		std::cout << "found at iter: " << stage_2_iter << std::endl;
		std::cout << "curr_normalized_integrality_gap: " << curr_normalized_integrality_gap_ << std::endl;
		std::cout << "cost: " << *(problem_->ComputeSolutionValue(curr_int_basis_)) << std::endl;
		for (int var_index = curr_int_basis_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = curr_int_basis_.find_next(var_index))
			std::cout << var_index << " ";

		std::cout << std::endl;
	}
	else
	{
		std::cout << "Não achou no limite de exec. Iter: " << stage_2_iter << std::endl;
		std::cout << "num active bin vars: " << problem_->curr_active_binary_vars().count() << std::endl;
	}

	(solution_)
		.time_stage2_ = timer->CurrentElapsedTime(ti);

	(solution_).num_iterations_stage2_ = stage_2_iter + (!K_SOLVE_STAGE_1 ? stage_0_iter : 0);
	(solution_).num_perturbations_stage1_ = num_perturbations_stage1;
	(solution_).num_perturbations_stage2_ = num_perturbations_stage2;
	(solution_).num_restarts_stage1_ = num_restarts_stage1;
	(solution_).num_restarts_stage2_ = num_restarts_stage2;

	if (found_int_basis_)
	{
		solution_.is_feasible_ = true;
		solution_.found_integer_ = true;
		solution_.bitset_vars_ = curr_int_basis_;
	}

	delete (ti);
	ti = nullptr;
	// std::cout << "iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
	// if(found_int_y_) std::cout << " * Y INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;
	// if(found_int_x_) std::cout << " * X INTEIRA!" << "   iterations: " << stage_0_iter + stage_1_iter + stage_2_iter << std::endl;

	// getchar(); getchar();
	return true;
}

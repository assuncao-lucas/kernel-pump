#pragma once

#include <ilcplex/ilocplex.h>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <utility>
#include <list>
#include "src/heuristic_solution.h"
#include "src/problem.h"

class FeasibilityPump
{
public:
	FeasibilityPump() = default;
	~FeasibilityPump();
	void Init(Problem *problem);
	bool Run(bool reset_fp_initial_basis_at_new_loop, double time_limit);
	boost::dynamic_bitset<> best_basis() const { return best_basis_; }
	double best_normalized_integrality_gap() const { return best_normalized_integrality_gap_; }
	bool found_int_basis() const { return found_int_basis_; }

private:
	double curr_alpha_ = 0.0; /* ratio by which the new obj is combined from distance and original obj */
	double previous_alpha_ = std::numeric_limits<double>::infinity();
	double previous_normalized_integrality_gap_ = std::numeric_limits<double>::infinity();
	double curr_normalized_integrality_gap_ = std::numeric_limits<double>::infinity();
	double best_normalized_integrality_gap_ = std::numeric_limits<double>::infinity();
	bool found_int_basis_ = false;

	boost::dynamic_bitset<> best_basis_;
	boost::dynamic_bitset<> curr_int_basis_;
	boost::dynamic_bitset<> previous_int_basis_;
	IloNumArray *curr_relax_basis_ = nullptr;

	std::list<std::pair<int, double>> curr_integrality_gaps_;

	Problem *problem_ = nullptr;

	void Reset();
	void RetrieveAndRoundBinaryVarsValues();
	void SetNewObjStage();

	FPHeuristicSolution solution_;
};

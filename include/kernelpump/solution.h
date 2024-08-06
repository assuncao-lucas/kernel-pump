#pragma once
#include <limits>

class Solution
{
public:
	explicit Solution() = default;
	virtual ~Solution() = default;

	bool is_feasible_ = false;
	double value_ = 0.0;
	double reopt_value_ = 0.0;
	double real_integrality_gap_ = std::numeric_limits<double>::infinity();		  // this value considers the original problem (not the presolved), on the best (possibly fraction) fractional solution obtained. It considers the distance to the nearest feasible integer value of each variable.
	double projection_integrality_gap_ = std::numeric_limits<double>::infinity(); // this value might not be reliable, cause it is based on the (possibly propagation-based) projection rouding of the (possibly presolved) problem.
	double total_time_spent_ = 0.0;
	int num_iterations_ = 0;
	double time_spent_building_kernel_buckets_ = 0.0;
	int num_buckets_ = 0;
	int last_bucket_visited_ = -1;
	int first_bucket_to_iter_pump_ = -1;
	int num_frac_ = 0; // number of integer/binary variables that are fractional in solution.
	int num_binary_vars_added_ = -1;
	int num_binary_vars_with_value_one_ = -1;

	void WriteToFile(std::string folder, std::string config_name, std::string instance_name, uint64_t seed) const;
};
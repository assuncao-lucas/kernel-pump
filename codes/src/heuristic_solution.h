#pragma once

#include <list>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "general.h"

class HeuristicSolution
{
public:
	explicit HeuristicSolution() = default;
	explicit HeuristicSolution(size_t num_vars);
	virtual ~HeuristicSolution();
	virtual void Reset(size_t num_vars);
	bool is_infeasible_ = false;
	bool is_feasible_ = false;
	bool is_optimal_ = false;
	double cost_ = 0.0;
	double total_time_spent_ = 0.0;

	bool operator==(HeuristicSolution &);

	friend std::ostream &operator<<(std::ostream &out, HeuristicSolution &sol);

	virtual void WriteToFile(std::string algo, std::string folder, std::string file_name) const;

	size_t num_vars_ = 0;
	boost::dynamic_bitset<> bitset_vars_;
};

class FPHeuristicSolution : public HeuristicSolution
{
public:
	explicit FPHeuristicSolution() = default;
	explicit FPHeuristicSolution(size_t num_vars);
	~FPHeuristicSolution() = default;
	int num_iterations_stage1_ = 0;
	int num_iterations_stage2_ = 0;
	int num_perturbations_stage1_ = 0;
	int num_perturbations_stage2_ = 0;
	int num_restarts_stage1_ = 0;
	int num_restarts_stage2_ = 0;
	bool found_integer_ = false;
	double time_stage1_ = 0;
	double time_stage2_ = 0;

	virtual void Reset(size_t num_vars);
	void WriteToFile(std::string algo, std::string folder, std::string file_name) const;
};

class KPHeuristicSolution : public HeuristicSolution
{
public:
	explicit KPHeuristicSolution() = default;
	explicit KPHeuristicSolution(int num_vars);
	~KPHeuristicSolution() = default;
	double time_spent_building_kernel_buckets_ = 0.0;
	bool found_integer_ = false;

	virtual void Reset(int num_vars);
	static std::string GenerateFileName(int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor);
	void WriteToFile(std::string algo, std::string folder, std::string file_name) const;
};
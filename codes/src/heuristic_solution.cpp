#include <fstream>
#include <sstream>
#include <iomanip>
#include <boost/dynamic_bitset.hpp>
#include "heuristic_solution.h"

HeuristicSolution::~HeuristicSolution()
{
}

HeuristicSolution::HeuristicSolution(size_t num_vars)
{
	HeuristicSolution::Reset(num_vars);
}

void HeuristicSolution::Reset(size_t num_vars)
{
	cost_ = 0;
	num_vars_ = num_vars;
	is_infeasible_ = false;
	is_feasible_ = false;
	is_optimal_ = false;
	total_time_spent_ = 0.0;

	bitset_vars_ = boost::dynamic_bitset<>(num_vars, 0);
}

bool HeuristicSolution::operator==(HeuristicSolution &other)
{
	return (bitset_vars_ == other.bitset_vars_);
}

std::ostream &operator<<(std::ostream &out, HeuristicSolution &sol)
{

	return out;
}

void HeuristicSolution::WriteToFile(std::string algo, std::string folder, std::string file_name) const
{
	// std::fstream file;
	// std::string path = "..//solutions//";
	// path.append(folder);
	// // struct stat sb;
	// // if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	// path.append("s_");
	// path.append(algo);
	// path.append("_");
	// path.append(file_name);
	// // std::cout << path << std::endl;

	// file.open(path.c_str(), std::fstream::out | std::fstream::app);

	// file << std::setprecision(5) << std::fixed;

	// file << K_FILE_DELIMITER << std::endl
	// 	 << "Total profit sum: " << profits_sum_ << std::endl
	// 	 << "Num routes: ";
	// if ((is_infeasible_) || (!(is_feasible_)))
	// 	file << "0" << std::endl;
	// else
	// {
	// 	file << num_routes_ << std::endl;

	// 	for (int i = 0; i < num_routes_; i++)
	// 	{
	// 		const Route &curr_route = (routes_vec_)[i];
	// 		file << curr_route.sum_profits_ << " " << curr_route.time_ << " ";

	// 		for (auto it = (curr_route.vertices_).begin(); it != (curr_route.vertices_).end(); ++it)
	// 		{
	// 			file << instance.getOriginalVertexPosition(*it) << " ";
	// 		}

	// 		// file << num_vertices_ - 1;
	// 		file << std::endl;
	// 	}
	// }

	// file.close();
}

FPHeuristicSolution::FPHeuristicSolution(size_t num_vars) : HeuristicSolution(num_vars)
{
}

void FPHeuristicSolution::Reset(size_t num_vars)
{
	HeuristicSolution::Reset(num_vars);

	num_iterations_stage1_ = 0;
	num_iterations_stage2_ = 0;
	num_perturbations_stage1_ = 0;
	num_perturbations_stage2_ = 0;
	num_restarts_stage1_ = 0;
	num_restarts_stage2_ = 0;
	found_integer_ = false;
	time_stage1_ = 0.0;
	time_stage2_ = 0.0;
}

void FPHeuristicSolution::WriteToFile(std::string algo, std::string folder, std::string file_name) const
{
	std::fstream file;
	std::string path = "..//solutions//";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	file.open(path.c_str(), std::fstream::out);

	file << std::setprecision(5) << std::fixed;

	if (is_infeasible_)
		file << "STATUS: INFEASIBLE" << std::endl;
	else if (found_integer_)
		file << "STATUS: FOUND INTEGER FEASIBLE" << std::endl;
	else
		file << "STATUS: FAILED TO FIND A FEASIBLE SOLUTION" << std::endl;
	file << "profit sum: " << cost_ << std::endl;
	file << "STAGE 1: " << std::endl;
	file << "# iterations: " << num_iterations_stage1_ << std::endl;
	file << "# perturbations: " << num_perturbations_stage1_ << std::endl;
	file << "# restarts: " << num_restarts_stage1_ << std::endl;
	file << "time(s): " << time_stage1_ << std::endl;
	file << "STAGE 2: " << std::endl;
	file << "# iterations: " << num_iterations_stage2_ << std::endl;
	file << "# perturbations: " << num_perturbations_stage2_ << std::endl;
	file << "# restarts: " << num_restarts_stage2_ << std::endl;
	file << "time(s): " << time_stage2_ << std::endl;

	file.close();

	HeuristicSolution::WriteToFile(algo, folder, file_name);
}

KPHeuristicSolution::KPHeuristicSolution(int num_vars) : HeuristicSolution(num_vars)
{
}

void KPHeuristicSolution::Reset(int num_vars)
{
	HeuristicSolution::Reset(num_vars);

	found_integer_ = false;
	time_spent_building_kernel_buckets_ = 0.0;
}

std::string KPHeuristicSolution::GenerateFileName(int ks_max_size_bucket, int ks_min_time_limit, int ks_max_time_limit, double ks_decay_factor)
{
	std::stringstream ss_decay_factor;
	ss_decay_factor << std::fixed << std::setprecision(2) << ks_decay_factor;
	std::string file_name = "kp_b" + std::to_string(ks_max_size_bucket) + "_[" + std::to_string(ks_max_time_limit) + "," + std::to_string(ks_min_time_limit) + "]_d" + ss_decay_factor.str();
	return file_name;
}

void KPHeuristicSolution::WriteToFile(std::string algo, std::string folder, std::string file_name) const
{
	std::fstream file;
	std::string path = "..//solutions//";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	file.open(path.c_str(), std::fstream::out);

	file << std::setprecision(5) << std::fixed;

	if (is_infeasible_)
		file << "STATUS: INFEASIBLE" << std::endl;
	else if (found_integer_)
		file << "STATUS: FOUND INTEGER FEASIBLE" << std::endl;
	else
		file << "STATUS: FAILED TO FIND A FEASIBLE SOLUTION" << std::endl;

	file << "time building kernel and buckets (s): " << time_spent_building_kernel_buckets_ << std::endl;
	file << "total time (s): " << total_time_spent_ << std::endl;

	file.close();

	HeuristicSolution::WriteToFile(algo, folder, file_name);
}
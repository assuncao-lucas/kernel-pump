#include <fstream>
#include <iomanip>
#include "kernelpump/solution.h"

void Solution::WriteToFile(std::string folder, std::string config_name, std::string instance_name, uint64_t seed) const
{
	std::fstream file;
	std::string path = "..//solutions//";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("//s_");
	path.append(config_name);
	path.append("_");
	path.append(instance_name);
	path.append("_");
	path.append(std::to_string(seed));
	path.append(".sol");
	// std::cout << path << std::endl;

	file.open(path.c_str(), std::fstream::out);

	file << std::setprecision(6) << std::fixed;

	if (is_feasible_)
		file << "STATUS: FOUND INTEGER FEASIBLE" << std::endl;
	else
		file << "STATUS: FAILED TO FIND AN INTEGER FEASIBLE SOLUTION" << std::endl;

	file << "time building kernel and buckets (s): " << time_spent_building_kernel_buckets_ << std::endl
		 << "total time (s): " << total_time_spent_ << std::endl
		 << "# iterations: " << num_iterations_ << std::endl
		 << "# buckets: " << num_buckets_ << std::endl
		 << "last bucket visited: " << last_bucket_visited_ << std::endl
		 << "first bucket to iter pump: " << first_bucket_to_iter_pump_ << std::endl
		 << "value: " << value_ << std::endl
		 << "reopt value: " << reopt_value_ << std::endl
		 << "real integrality gap: " << real_integrality_gap_ << std::endl
		 << "projection integrality gap: " << projection_integrality_gap_ << std::endl
		 << "num frac: " << num_frac_;

	if (num_binary_vars_added_ != -1)
	{
		file << std::endl
			 << "num bin vars added: " << num_binary_vars_added_ << std::endl
			 << "num bin vars with value 1: " << num_binary_vars_with_value_one_ << std::endl;
	}

	file.close();
}
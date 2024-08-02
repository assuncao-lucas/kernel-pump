#include <iostream>
#include <algorithm>
#include <queue>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <dirent.h>
#include <unistd.h>
#include <getopt.h>
#include <cstdlib>
#include <boost/dynamic_bitset.hpp>
#include "src/graph.h"
#include "src/instance.h"
#include "src/exact/formulations.h"
#include "src/timer.h"
#include "src/general.h"
#include "src/graph_algorithms.h"
#include "src/solution.hpp"
#include "src/heuristic_solution.h"
#include "src/feasibility_pump/feasibility_pump.h"
#include "src/kernel_search/kernel_search.h"
#include "src/ALNS/ALNS.h"
#include "src/simulated_annealing/simulated_annealing.h"

// Instance* GenerateFilesFromOPinstancesIter(std::string dir_path, std::string file_name, double mandatory_percentage)
// {
// 	std::string curr_file = dir_path;
// 	curr_file.append(file_name);
// 	Graph * graph = NULL;
// 	Instance * new_instance = NULL;

// 	std::fstream file;
// 	file.open(curr_file.c_str(), std::fstream::in);
// 	if(!file.is_open()){
// 		std::cout << "Could not open file" << std::endl;
// 		throw 1;
// 		return NULL;
// 	}

// 	int num_vertices = 0, num_vehicles = 0, num_mandatory = 0;
// 	double limit = 0.0;

// 	/*std::cout << file_name << std::endl;
// 	std::cout << num_vertices << std::endl;
// 	std::cout << num_mandatory << std::endl;
// 	std::cout << num_vehicles << std::endl;
// 	std::cout << limit << std::endl;
// 	std::cout << "*******************" << std::endl;*/

// 	std::vector<std::pair<double,double>> coordinates;
// 	std::pair<double,double> curr_coordinate;
// 	int curr_profit;

// 	std::vector<int> * profits_vec = new std::vector<int>();

// 	file >> limit >> num_vehicles;

// 	// skips position 1 in coordinates vector: will be filled with last node (will be a mandatory node)
// 	while(file >> curr_coordinate.first)
// 	{
// 		num_vertices++;
// 		file >> curr_coordinate.second >> curr_profit;
// 		coordinates.push_back(curr_coordinate);
// 		profits_vec->push_back(curr_profit);
// 	}

// 	file.close();

// 	// relocates destination to the end of the vector
// 	curr_coordinate.first = coordinates[1].first;
// 	curr_coordinate.second = coordinates[1].second;
// 	curr_profit = (*profits_vec)[1];

// 	coordinates[1].first = coordinates[num_vertices-1].first;
// 	coordinates[1].second = coordinates[num_vertices-1].second;
// 	(*profits_vec)[1] = (*profits_vec)[num_vertices-1];

// 	coordinates[num_vertices-1].first = curr_coordinate.first;
// 	coordinates[num_vertices-1].second = curr_coordinate.second;
// 	(*profits_vec)[num_vertices-1] = curr_profit;

// 	num_mandatory = (int)ceil(mandatory_percentage * num_vertices);

// 	graph = new Graph(num_vertices,num_mandatory,(&((*profits_vec)[0])));

// 	/*for(int i = 0; i < num_vertices; i++)
// 	{
// 	std::cout << i << " : (" << coordinates[i].first << ", " << coordinates[i].second << ")";
// 	if(i <= num_mandatory) std::cout << " | profit: 0.0" << std::endl;
// 	else std::cout << " | profit: " << profits[i - num_mandatory-1] << std::endl;
// 	getchar();
// }*/

// for(int i = 0; i < num_vertices; i++)
// {
// 	for(int j = i+1; j < num_vertices; j ++)
// 	{
// 		graph->AddEdge(i,j,euclidian_distance(coordinates[i], coordinates[j]));
// 	}
// }

// new_instance = new Instance(graph, num_vehicles,limit);

// std::vector<bool> selected_vertices(num_vertices, false);
// int iter_mandatory = 0, cont_mandatory = 0;
// while(cont_mandatory < num_mandatory)
// {
// 	iter_mandatory = rand()%(num_vertices-2) + 1;
// 	if(!selected_vertices[iter_mandatory])
// 	{
// 		selected_vertices[iter_mandatory] = true;
// 		cont_mandatory++;
// 		new_instance->mandatory_list_.push_back(iter_mandatory);
// 	}
// }

// return new_instance;
// }

// Instance* GenerateFilesFromTOPinstancesIter(std::string dir_path, std::string file_name, double mandatory_percentage)
// {
// 	std::string curr_file = dir_path;
// 	curr_file.append(file_name);
// 	Graph * graph = NULL;
// 	Instance * new_instance = NULL;

// 	std::cout << file_name << std::endl;

// 	std::fstream file;
// 	file.open(curr_file.c_str(), std::fstream::in);
// 	if(!file.is_open()){
// 		std::cout << "Could not open file" << std::endl;
// 		throw 1;
// 		return NULL;
// 	}

// 	int num_vertices = 0, num_vehicles = 0, num_mandatory = 0;
// 	double limit = 0.0;

// 	std::string line;
// 	getline(file,line);

// 	std::size_t pos_1 = line.find_first_of(" ");

// 	std::stringstream parameter;
// 	parameter << line.substr(pos_1 + 1);
// 	parameter >> num_vertices;

// 	num_mandatory = (int)ceil(mandatory_percentage * num_vertices);

// 	getline(file,line);
// 	pos_1 = line.find_first_of(" ");

// 	parameter.clear();
// 	parameter << line.substr(pos_1 + 1);
// 	parameter >> num_vehicles;

// 	getline(file,line);
// 	pos_1 = line.find_first_of(" ");

// 	parameter.clear();
// 	parameter << line.substr(pos_1 + 1);
// 	parameter >> limit;

// 	/*std::cout << file_name << std::endl;
// 	std::cout << num_vertices << std::endl;
// 	std::cout << num_mandatory << std::endl;
// 	std::cout << num_vehicles << std::endl;
// 	std::cout << limit << std::endl;
// 	std::cout << "*******************" << std::endl;*/

// 	std::vector<std::pair<double,double>> coordinates(num_vertices);
// 	int * profits = new int[num_vertices];

// 	// skips position 1 in coordinates vector: will be filled with last node (will be a mandatory node)
// 	for(int i = 0; i < num_vertices; i++)
// 	{
// 		file >> (coordinates[i]).first >> (coordinates[i]).second >> profits[i];
// 		//std::cout << (coordinates[i]).first << " " << (coordinates[i]).second << " " << profits[i];
// 		//getchar();getchar();
// 	}

// 	file.close();

// 	graph = new Graph(num_vertices,num_mandatory,profits);

// 	/*for(int i = 0; i < num_vertices; i++)
// 	{
// 	std::cout << i << " : (" << coordinates[i].first << ", " << coordinates[i].second << ")";
// 	if(i <= num_mandatory) std::cout << " | profit: 0.0" << std::endl;
// 	else std::cout << " | profit: " << profits[i - num_mandatory-1] << std::endl;
// 	getchar();
// }*/

// for(int i = 0; i < num_vertices; i++)
// {
// 	for(int j = i+1; j < num_vertices; j ++)
// 	{
// 		graph->AddEdge(i,j,euclidian_distance(coordinates[i], coordinates[j]));
// 	}
// }

// new_instance = new Instance(graph, num_vehicles,limit);

// std::vector<bool> selected_vertices(num_vertices, false);
// int iter_mandatory = 0, cont_mandatory = 0;
// while(cont_mandatory < num_mandatory)
// {
// 	iter_mandatory = rand()%(num_vertices-2) + 1;
// 	if(!selected_vertices[iter_mandatory])
// 	{
// 		selected_vertices[iter_mandatory] = true;
// 		cont_mandatory++;
// 		new_instance->mandatory_list_.push_back(iter_mandatory);
// 	}
// }

// return new_instance;
// }

// int GenerateFilesFromOPinstances(std::string dir_path, double mandatory_percentage)
// {
// 	Instance * inst = NULL;
// 	DIR *dir;
// 	struct dirent *ent;
// 	if ((dir = opendir (dir_path.c_str())) != NULL)
// 	{
// 		/* print all the files and directories within directory */
// 		while ((ent = readdir (dir)) != NULL)
// 		{
// 			std::string curr_file(ent->d_name);
// 			if((curr_file != ".")&&(curr_file != "..")&&( (curr_file.size() > 0)&& (curr_file[curr_file.size() -1] != '~') )) inst = GenerateFilesFromOPinstancesIter(dir_path, curr_file, mandatory_percentage);
// 			if(inst != NULL)
// 			{
// 				std::string folder;
// 				std::stringstream s_percentage;
// 				std::string percentage;
// 				std::fstream output;

// 				std::size_t pos = dir_path.find_first_of("OP");

// 				folder = dir_path.substr(0,pos);
// 				folder.append("ST");
// 				folder.append(dir_path.substr(pos));
// 				folder.append("Set_");

// 				s_percentage << mandatory_percentage;
// 				percentage = s_percentage.str();

// 				folder.append(percentage);
// 				folder.append("//");
// 				folder.append(curr_file);
// 				std::cout << folder << std::endl;
// 				inst->WriteToFile(folder,curr_file,mandatory_percentage);
// 				delete inst;
// 				inst = NULL;
// 			}
// 		}
// 		closedir (dir);
// 	}else
// 	{
// 		return -1;
// 	}
// 	return 0;
// }

// int GenerateFilesFromTOPinstances(std::string dir_path, double mandatory_percentage, bool repair = false)
// {
// 	Instance * inst = NULL;
// 	DIR *dir;

// 	struct dirent *ent;
// 	if ((dir = opendir (dir_path.c_str())) != NULL)
// 	{
// 		/* print all the files and directories within directory */
// 		while ((ent = readdir (dir)) != NULL)
// 		{
// 			std::string curr_file(ent->d_name);
// 			if((curr_file != ".")&&(curr_file != "..")&&( (curr_file.size() > 0)&& (curr_file[curr_file.size() -1] != '~') )) inst = GenerateFilesFromTOPinstancesIter(dir_path, curr_file, mandatory_percentage);
// 			if(inst != NULL)
// 			{
// 				std::string folder;
// 				std::stringstream s_percentage;
// 				std::string percentage;
// 				std::fstream output;

// 				std::size_t pos = dir_path.find_first_of("TOP");

// 				folder = dir_path.substr(0,pos);
// 				folder.append("S");
// 				folder.append(dir_path.substr(pos));
// 				folder.append("Set_");

// 				s_percentage << mandatory_percentage;
// 				percentage = s_percentage.str();

// 				folder.append(percentage);
// 				folder.append("//");
// 				folder.append(curr_file);

// 				std::string folder2 = folder;
// 				folder2.replace(14,4,"STOP_backup");

// 				//std::cout << folder << std::endl;
// 				//std::cout << folder2 << std::endl;

// 				if(repair)
// 				{
// 					Instance inst_old(folder2);

// 					inst->mandatory_list_ = inst_old.mandatory_list_;
// 				}

// 				inst->WriteToFile(folder,curr_file,mandatory_percentage);
// 				//getchar(); getchar();
// 				delete inst;
// 				inst = NULL;
// 			}
// 		}
// 		closedir (dir);
// 	}else
// 	{
// 		return -1;
// 	}
// 	return 0;
// }

void GenrateRelaxCSVTable(std::vector<std::string> dirs, bool stop, double total_time_limit)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_none");
	/*algorithms.push_back("relax_all_gccs");
	//algorithms.push_back("relax_gccs_cb_csc3");
	algorithms.push_back("relax_all_ccs");
	//algorithms.push_back("relax_ccs_cb_csc3");
	algorithms.push_back("relax_all_gccs_ccs");
	//algorithms.push_back("relax_gccs_ccs_cb_csc3");
	//algorithms.push_back("relax_all_gccs_ccs_lcis");
	//algorithms.push_back("relax_clique_ccs_active_vertices");
	algorithms.push_back("relax_clique_ccs_active_vertices_unilateral");
	//algorithms.push_back("relax_clique_ccs_total_active_vertices_unilateral");
	algorithms.push_back("relax_maximum_clique_ccs_active_vertices_unilateral");
	//algorithms.push_back("relax_lcis_clique_ccs_active_vertices_unilateral");

	algorithms.push_back("relax_angle_0.1_maximum_clique_ccs_active_vertices_unilateral");
	algorithms.push_back("relax_angle_0.05_clique_ccs_active_vertices_unilateral");*/
	algorithms.push_back("cb_csc3");
	algorithms.push_back("relax_angle_0.03_clique_ccs_active_vertices_unilateral");

	std::fstream output;
	std::string output_name;
	if (stop)
		output_name = ".//tables//CSV//stop_relax_table.csc";
	else
		output_name = ".//tables//CSV//relax_table.csv";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_time_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_cuts_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_time(algorithms.size(), 0.0);
	std::vector<double> total_avg_cuts(algorithms.size(), 0.0);
	int total_num_instances = 0;

	output << "set";
	for (size_t j = 0; j < algorithms.size(); j++)
	{
		output << ";" << algorithms[j] << ";" << algorithms[j];
	}
	output << std::endl;
	for (size_t j = 0; j < algorithms.size(); j++)
	{
		output << ";tempo(s);#cuts";
	}
	output << std::endl;
	for (size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> time_per_algo(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> cuts_per_algo(algorithms.size(), std::vector<double>());
		std::vector<double> avg_time(algorithms.size(), 0.0);
		std::vector<double> avg_cuts(algorithms.size(), 0.0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);

		total_num_instances += (instances.size());
		for (size_t i = 0; i < instances.size(); i++)
		{
			for (size_t j = 0; j < algorithms.size(); j++)
			{
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_time;
				std::string status;
				double time = 0.0;
				int num_cuts = 0;
				std::string line;

				for (int i = 0; i < 5; i++)
					getline(input, line);

				size_t pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);
				s_time >> time;
				// std::cout << line << "    " << time << std::endl;
				// getchar(); getchar();

				for (int i = 0; i < 4; i++)
					getline(input, line);

				for (int i = 0; i < K_NUM_TYPES_CALLBACKS; i++)
				{
					std::stringstream s_cuts;
					int cuts_iter = 0;
					getline(input, line);
					pos = line.find_first_of(":");
					size_t pos2 = line.find_first_of("/");
					s_cuts << line.substr(pos + 2, pos2 - pos - 2);
					s_cuts >> cuts_iter;
					num_cuts += cuts_iter;

					// std::cout << line << "    " << cuts_iter << std::endl;
					// getchar(); getchar();
				}
				time_per_algo[j].push_back(time);
				total_time_per_algo[j].push_back(time);
				total_avg_time[j] += time;
				avg_time[j] += time;

				// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				cuts_per_algo[j].push_back(num_cuts);
				total_cuts_per_algo[j].push_back(num_cuts);

				total_avg_cuts[j] += num_cuts;

				avg_cuts[j] += num_cuts;
				input.close();
			}

			// getchar();getchar();
		}

		output << j + 1;

		for (size_t j = 0; j < algorithms.size(); j++)
		{
			if ((time_per_algo[j]).size() > 0)
				avg_time[j] /= (1.0 * ((time_per_algo[j]).size()));
			else
				avg_time[j] = -1;
			avg_cuts[j] /= (1.0 * ((cuts_per_algo[j]).size()));
			output << ";" << avg_time[j] << ";" << avg_cuts[j];
		}
		output << std::endl;
	}

	output << "Total";
	for (size_t j = 0; j < algorithms.size(); j++)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_time[j] /= (1.0 * ((total_time_per_algo[j]).size()));
		total_avg_cuts[j] /= (1.0 * ((total_cuts_per_algo[j]).size()));
		output << ";" << total_avg_time[j] << ";" << total_avg_cuts[j];
	}

	output << std::endl;

	output.close();
}

void GenerateAlgorithmsLatexTable(std::vector<std::string> dirs, bool stop, double total_time_limit)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	// algorithms.push_back("bc_b1");
	// algorithms.push_back("csc3");
	// algorithms.push_back("csc3_2");
	// algorithms.push_back("cb_csc3");
	// algorithms.push_back("cb2_csc3");

	/*std::string algo = "cb_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	algo = "cb_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);*/

	std::string algo = "cb_csc3_";
	if (stop)
		algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	algo = "cb2_csc3_";
	if (stop)
		algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	std::fstream output;
	std::string output_name;
	if (stop)
		output_name = ".//tables//latex//stop_table_algorithms_new.txt";
	else
		output_name = ".//tables//latex//table_algorithms_new.txt";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_time_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_time(algorithms.size(), 0.0);
	std::vector<double> total_avg_gap(algorithms.size(), 0.0);
	std::vector<int> total_num_optimal(algorithms.size(), 0);
	int total_num_instances = 0;

	for (size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> time_per_algo(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo(algorithms.size(), std::vector<double>());
		std::vector<double> avg_time(algorithms.size(), 0.0);
		std::vector<double> avg_gap(algorithms.size(), 0.0);
		std::vector<double> st_dev(algorithms.size(), 0.0);

		std::vector<int> num_optimal(algorithms.size(), 0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);

		total_num_instances += (instances.size());
		for (size_t i = 0; i < instances.size(); i++)
		{
			for (size_t j = 0; j < algorithms.size(); j++)
			{
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lb, s_ub, s_time;
				std::string status;
				double lb = 0.0, ub = 0.0, time = 0.0, gap = 0.0;
				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_lb << line.substr(pos + 2);
				if (s_lb.str() == "-inf")
					lb = -1;
				else
					s_lb >> lb;

				getline(input, line);
				pos = line.find_first_of(":");
				s_ub << line.substr(pos + 2);
				if (s_ub.str() == "inf")
					ub = -1;
				else
					s_ub >> ub;

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);
				s_time >> time;

				if ((status != "OPTIMAL") && (status != "INFEASIBLE"))
				{
					if ((!(double_equals(lb, -1))) && (!(double_equals(ub, -1))))
					{
						if (double_less(time, total_time_limit))
						{
							(num_optimal[j])++;
							(total_num_optimal[j])++;
							time_per_algo[j].push_back(time);
							total_time_per_algo[j].push_back(time);

							total_avg_time[j] += time;
							avg_time[j] += time;
						}
						else
							gap = (100.0 * (ub - lb)) / ub;
					}
					else
						gap = 100.0;
				}
				else
				{
					(num_optimal[j])++;
					(total_num_optimal[j])++;
					time_per_algo[j].push_back(time);
					total_time_per_algo[j].push_back(time);

					total_avg_time[j] += time;
					avg_time[j] += time;
				}

				// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				if (!double_equals(gap, 0.0))
				{
					gap_per_algo[j].push_back(gap);
					total_gap_per_algo[j].push_back(gap);

					total_avg_gap[j] += gap;

					avg_gap[j] += gap;
				}
				input.close();
			}

			// getchar();getchar();
		}

		output << j + 1;

		if (stop)
			output << "\\_5\\%";

		for (size_t j = 0; j < algorithms.size(); j++)
		{
			if ((time_per_algo[j]).size() > 0)
				avg_time[j] /= (1.0 * ((time_per_algo[j]).size()));
			else
				avg_time[j] = -1;
			if ((gap_per_algo[j]).size() > 0)
			{
				avg_gap[j] /= (1.0 * ((gap_per_algo[j]).size()));
				st_dev[j] = StDev(gap_per_algo[j], avg_gap[j]);
			}
			else
				avg_gap[j] = st_dev[j] = -1;
			output << " & & " << num_optimal[j] << "/" << instances.size() << " & " << avg_time[j] << " & " << avg_gap[j] << " & " << st_dev[j];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for (size_t j = 0; j < algorithms.size(); j++)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		if ((total_time_per_algo[j]).size() > 0)
			total_avg_time[j] /= (1.0 * ((total_time_per_algo[j]).size()));
		else
			total_avg_time[j] = -1;
		if ((total_gap_per_algo[j]).size() > 0)
			total_avg_gap[j] /= (1.0 * ((total_gap_per_algo[j]).size()));
		else
			total_avg_gap[j] = -1;
		output << "& & " << total_num_optimal[j] << "/" << total_num_instances << " & " << total_avg_time[j] << " & " << total_avg_gap[j] << " & " << StDev(total_gap_per_algo[j], total_avg_gap[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateValidIneqsActiveLatexTable(std::vector<std::string> dirs, bool stop)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_b1");
	algorithms.push_back("relax_none_cb_csc3");

	std::fstream output;
	std::string output_name;
	if (stop)
		output_name = ".//tables//latex//stop_table_valid_ineq_active.txt";
	else
		output_name = ".//tables//latex//table_valid_ineq_active.txt";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_percentage_active_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_percentage_active(algorithms.size(), 0.0);

	std::vector<int> total_num_active_per_algo(algorithms.size(), 0);
	std::vector<int> total_count(algorithms.size(), 0);

	for (size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> percentage_active_per_algo(algorithms.size(), std::vector<double>());
		std::vector<double> avg_percentage_active(algorithms.size(), 0.0);
		std::vector<double> st_dev(algorithms.size(), 0.0);

		std::vector<int> num_active_per_algo(algorithms.size(), 0);
		std::vector<int> count(algorithms.size(), 0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);
		for (size_t i = 0; i < instances.size(); i++)
		{
			// std::cout << "* " << instances[i] << std::endl;
			int num_active = 0, num_total = 0;
			double percentage_active = 0.0;
			for (size_t j = 0; j < algorithms.size(); j++)
			{
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_active, s_total;
				int curr_active = 0, curr_total = 0;
				std::string line;

				for (int k = 0; k < 12; k++)
					getline(input, line);
				size_t pos = line.find_first_of(":"), pos2 = line.find_first_of("/");

				s_active << line.substr(pos + 2, pos2 - pos - 2);
				s_active >> num_active;

				s_total << line.substr(pos2 + 1);
				s_total >> num_total;

				// std::cout << num_active << "/" << num_total << std::endl;
				// getchar();getchar();

				if (num_total != 0)
				{
					percentage_active = (num_active * 1.0) / num_total;
					percentage_active_per_algo[j].push_back(percentage_active);
					total_percentage_active_per_algo[j].push_back(percentage_active);
					total_avg_percentage_active[j] += percentage_active;
					avg_percentage_active[j] += percentage_active;

					(num_active_per_algo[j]) += num_active;
					(total_num_active_per_algo[j]) += num_active;
					(count[j])++;
					(total_count[j])++;
				}

				// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;
				input.close();
			}
		}

		output << j + 1;

		for (size_t j = 0; j < algorithms.size(); j++)
		{
			avg_percentage_active[j] /= (1.0 * ((percentage_active_per_algo[j]).size()));
			st_dev[j] = StDev(percentage_active_per_algo[j], avg_percentage_active[j]);
			output << " & & " << (1.0 * num_active_per_algo[j]) / (count[j]) << " & " << avg_percentage_active[j] << " & " << st_dev[j];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";

	for (size_t j = 0; j < algorithms.size(); j++)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_percentage_active[j] /= (1.0 * ((total_percentage_active_per_algo[j]).size()));
		output << "& & " << (1.0 * total_num_active_per_algo[j]) / (total_count[j]) << " & " << total_avg_percentage_active[j] << " & " << StDev(total_percentage_active_per_algo[j], total_avg_percentage_active[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateValidIneqsImprovementsLatexTable(std::vector<std::string> dirs, bool stop)
{
	std::string curr_file;
	std::vector<std::string> algorithms;
	std::vector<std::string> algorithms2;

	algorithms.push_back("relax_cb_csc3");
	algorithms.push_back("relax_FBs_cb_csc3");

	algorithms2.push_back("relax_b1_sem");
	algorithms2.push_back("relax_b1");

	std::fstream output;
	std::string output_name;
	if (stop)
		output_name = ".//tables//latex//stop_table_valid_ineq_improvements.txt";
	else
		output_name = ".//tables//latex//table_valid_ineq_improvements.txt";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_improvement_per_algo(algorithms.size(), std::vector<double>()), total_improvement_per_algo2(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_improvement(algorithms.size(), 0.0), total_avg_improvement2(algorithms.size(), 0.0);

	for (size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> improvement_per_algo(algorithms.size(), std::vector<double>()), improvement_per_algo2(algorithms.size(), std::vector<double>());
		std::vector<double> avg_improvement(algorithms.size(), 0.0), avg_improvement2(algorithms.size(), 0.0);
		std::vector<double> st_dev(algorithms.size(), 0.0), st_dev2(algorithms.size(), 0.0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);
		for (size_t i = 0; i < instances.size(); i++)
		{
			double original_lp = 0.0, original_lp2 = 0.0;
			;
			for (size_t j = 0; j < algorithms.size(); j++)
			{
				double curr_improvement = 0.0;
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp;
				std::string status;
				double lp = 0.0;
				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);
				s_lp >> lp;

				if (j == 0)
				{
					if ((double_equals(lp, -1)) || (status == "INFEASIBLE"))
					{
						original_lp = -1.0;
					}
					else
						original_lp = lp;
					curr_improvement = 0.0;
				}
				else
				{
					if ((double_equals(lp, -1)) || (status == "INFEASIBLE") || (double_equals(original_lp, -1)))
					{
						curr_improvement = 0.0;
					}
					else
					{
						if (double_equals(original_lp, 0.0))
							curr_improvement = 0;
						else
							curr_improvement = (100 * (original_lp - lp)) / original_lp;
					}
				}

				// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				improvement_per_algo[j].push_back(curr_improvement);
				total_improvement_per_algo[j].push_back(curr_improvement);
				total_avg_improvement[j] += curr_improvement;
				avg_improvement[j] += curr_improvement;
				input.close();
			}

			for (size_t j = 0; j < algorithms2.size(); j++)
			{
				double curr_improvement = 0.0;
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms2[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp;
				std::string status;
				double lp = 0.0;
				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);
				s_lp >> lp;

				if (j == 0)
				{
					if ((double_equals(lp, -1)) || (status == "INFEASIBLE"))
					{
						original_lp2 = -1.0;
					}
					else
						original_lp2 = lp;
					curr_improvement = 0.0;
				}
				else
				{
					if ((double_equals(lp, -1)) || (status == "INFEASIBLE") || (double_equals(original_lp2, -1)))
					{
						curr_improvement = 0.0;
					}
					else
					{
						if (double_equals(original_lp2, 0.0))
							curr_improvement = 0;
						else
							curr_improvement = (100 * (original_lp2 - lp)) / original_lp2;
					}
				}

				// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				improvement_per_algo2[j].push_back(curr_improvement);
				total_improvement_per_algo2[j].push_back(curr_improvement);
				total_avg_improvement2[j] += curr_improvement;
				avg_improvement2[j] += curr_improvement;
				input.close();
			}

			// getchar();getchar();
		}

		output << j + 1;

		for (size_t j = 1; j < algorithms2.size(); j++)
		{
			avg_improvement2[j] /= (1.0 * (instances.size()));
			st_dev2[j] = StDev(improvement_per_algo2[j], avg_improvement2[j]);
			output << " & & " << avg_improvement2[j] << " & " << st_dev2[j];
		}
		for (size_t j = 1; j < algorithms.size(); j++)
		{
			avg_improvement[j] /= (1.0 * (instances.size()));
			st_dev[j] = StDev(improvement_per_algo[j], avg_improvement[j]);
			output << " & & " << avg_improvement[j] << " & " << st_dev[j];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for (size_t j = 1; j < algorithms2.size(); j++)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_improvement2[j] /= (1.0 * ((total_improvement_per_algo2[j]).size()));
		output << "& & " << total_avg_improvement2[j] << " & " << StDev(total_improvement_per_algo2[j], total_avg_improvement2[j]);
	}

	for (size_t j = 1; j < algorithms.size(); j++)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_improvement[j] /= (1.0 * ((total_improvement_per_algo[j]).size()));
		output << "& & " << total_avg_improvement[j] << " & " << StDev(total_improvement_per_algo[j], total_avg_improvement[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateGapsNodesCutsLatexTable(std::vector<std::string> dirs, bool stop, double time_limit)
{
	// std::cout << "3" << std::endl; getchar(); getchar();
	std::vector<std::string> algorithms;

	algorithms.push_back("bc_b1");
	algorithms.push_back("cb_csc3");
	algorithms.push_back("cb2_csc3");

	std::fstream output;
	std::string output_name;
	if (stop)
		output_name = ".//tables//latex//stop_table_gaps_nodes_cuts.txt";
	else
		output_name = ".//tables//latex//table_gaps_nodes_cuts.txt";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	for (size_t k = 0; k < algorithms.size(); k++)
	{
		std::vector<double> total_gaps;
		double total_avg_gaps = 0.0;
		std::vector<double> total_nodes;
		double total_avg_nodes = 0.0;
		std::vector<std::vector<double>> total_cuts(K_NUM_TYPES_CALLBACKS, std::vector<double>());
		std::vector<double> total_avg_cuts(K_NUM_TYPES_CALLBACKS, 0.0);

		output << algorithms[k] << std::endl;

		for (size_t j = 0; j < dirs.size(); j++)
		{
			std::vector<double> gaps;
			double avg_gaps = 0.0;
			std::vector<double> nodes;
			double avg_nodes = 0.0;
			std::vector<std::vector<double>> cuts(K_NUM_TYPES_CALLBACKS, std::vector<double>());
			std::vector<double> avg_cuts(K_NUM_TYPES_CALLBACKS, 0.0);

			std::vector<std::string> instances;
			std::string folder = dirs[j].substr(20);
			AddInstancesFromDirectory(dirs[j], instances, false);

			for (size_t i = 0; i < instances.size(); i++)
			{
				double lb = 0.0, ub = 0.0;
				double curr_gap = 0.0, curr_nodes = 0.0, curr_cuts = 0.0, time = 0.0;

				std::string curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[k]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;

				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_nodes, s_lb, s_ub, s_time;
				std::string status;

				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_lb << line.substr(pos + 2);
				if (s_lb.str() == "-inf")
					lb = -1;
				else
					s_lb >> lb;

				getline(input, line);
				pos = line.find_first_of(":");
				s_ub << line.substr(pos + 2);
				if (s_ub.str() == "inf")
					ub = -1;
				else
					s_ub >> ub;

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);
				s_time >> time;

				if ((status != "OPTIMAL") && (status != "INFEASIBLE"))
				{
					if ((!(double_equals(lb, -1))) && (!(double_equals(ub, -1))))
					{
						if (double_less(time, time_limit))
						{
							curr_gap = 0.0;
						}
						else
							curr_gap = (100.0 * (ub - lb)) / ub;
					}
					else
						curr_gap = 100.0;
				}
				else
				{
					curr_gap = 0.0;
				}

				getline(input, line);
				getline(input, line);
				getline(input, line);

				for (int l = 0; l < K_NUM_TYPES_CALLBACKS; ++l)
					getline(input, line);

				getline(input, line);
				getline(input, line);

				// read cuts here
				for (int l = 0; l < K_NUM_TYPES_CALLBACKS; ++l)
				{
					std::stringstream s_cuts;
					getline(input, line);
					pos = line.find_first_of(":");
					size_t pos2 = line.find_first_of("/");
					s_cuts << line.substr(pos + 2, pos2 - pos - 2);
					curr_cuts = 0.0;
					s_cuts >> curr_cuts;

					cuts[l].push_back(curr_cuts);
					total_cuts[l].push_back(curr_cuts);
					avg_cuts[l] += curr_cuts;
					total_avg_cuts[l] += curr_cuts;
					// std::cout << line << " " << curr_cuts << std::endl;
					// getchar();getchar();
				}

				getline(input, line);
				pos = line.find_first_of(":");
				s_nodes << line.substr(pos + 2);
				s_nodes >> curr_nodes;

				gaps.push_back(curr_gap);
				total_gaps.push_back(curr_gap);
				avg_gaps += curr_gap;
				total_avg_gaps += curr_gap;

				nodes.push_back(curr_nodes);
				total_nodes.push_back(curr_nodes);
				avg_nodes += curr_nodes;
				total_avg_nodes += curr_nodes;

				input.close();
			}
			// getchar();getchar();

			avg_gaps /= (gaps.size());
			avg_nodes /= (nodes.size());
			output << j + 1;
			if (stop)
				output << "\\_5\\%";

			output << " & & " << avg_gaps << " & " << StDev(gaps, avg_gaps) << " & & " << avg_nodes << " &";

			for (int l = 0; l < K_NUM_TYPES_CALLBACKS; ++l)
			{
				if (algorithms[k] == "bc_b1")
				{
					if (l == K_TYPE_GCC_CUT)
					{
						avg_cuts[l] /= (cuts[l].size());
						output << " & " << avg_cuts[l];
					}
				}

				if (algorithms[k] == "cb_csc3")
				{
					if ((l == K_TYPE_GCC_CUT) || (l == K_TYPE_CONFLICT_CUT) || (l == K_TYPE_COVER_BOUND_CUT))
					{
						avg_cuts[l] /= (cuts[l].size());
						output << " & " << avg_cuts[l];
					}
				}

				if (algorithms[k] == "cb2_csc3")
				{
					if ((l == K_TYPE_CLIQUE_CONFLICT_CUT) || (l == K_TYPE_COVER_BOUND_CUT))
					{
						avg_cuts[l] /= (cuts[l].size());
						output << " & " << avg_cuts[l];
					}
				}
			}

			output << "\\\\" << std::endl;
		}

		total_avg_gaps /= (total_gaps.size());
		total_avg_nodes /= (total_nodes.size());
		output << " & & " << total_avg_gaps << " & " << StDev(total_gaps, total_avg_gaps) << " & & " << total_avg_nodes << " &";

		for (int l = 0; l < K_NUM_TYPES_CALLBACKS; ++l)
		{
			if (algorithms[k] == "bc_b1")
			{
				if (l == K_TYPE_GCC_CUT)
				{
					total_avg_cuts[l] /= (total_cuts[l].size());
					output << " & " << total_avg_cuts[l];
				}
			}

			if (algorithms[k] == "cb_csc3")
			{
				if ((l == K_TYPE_GCC_CUT) || (l == K_TYPE_CONFLICT_CUT) || (l == K_TYPE_COVER_BOUND_CUT))
				{
					total_avg_cuts[l] /= (total_cuts[l].size());
					output << " & " << total_avg_cuts[l];
				}
			}

			if (algorithms[k] == "cb2_csc3")
			{
				if ((l == K_TYPE_CLIQUE_CONFLICT_CUT) || (l == K_TYPE_COVER_BOUND_CUT))
				{
					total_avg_cuts[l] /= (total_cuts[l].size());
					output << " & " << total_avg_cuts[l];
				}
			}
		}

		output << "\\\\" << std::endl;
	}
	output.close();
}

void GenerateRootPrimalAndDualGapsLatexTable(std::vector<std::string> dirs, bool stop, double time_limit, bool user_cuts)
{
	// std::cout << "3" << std::endl; getchar(); getchar();
	std::vector<std::pair<std::string, std::string>> algorithms;

	// algorithms.push_back(std::pair<std::string,std::string>("bc_b1_initial_primal_bound","bc_b1"));
	algorithms.push_back(std::pair<std::string, std::string>("cb_csc3_initial_primal_bound", "cb_csc3"));
	algorithms.push_back(std::pair<std::string, std::string>("cb2_csc3_initial_primal_bound", "cb2_csc3"));

	std::vector<std::pair<std::string, std::string>> algorithms2;
	// algorithms2.push_back(std::pair<std::string,std::string>("bc_b1_initial_primal_bound","bc_b1"));
	// algorithms2.push_back(std::pair<std::string,std::string>("cb_csc3_initial_primal_bound","cb_csc3"));
	algorithms2.push_back(std::pair<std::string, std::string>("cb2_csc3_initial_primal_bound", "cb2_csc3"));
	;

	std::fstream output;
	std::string output_name;
	if (stop)
		output_name = ".//tables//latex//stop_table_primal_dual_gaps";
	else
		output_name = ".//tables//latex//table_primal_dual_gaps";

	if (!user_cuts)
		output_name += "_no_user_cuts";

	output_name += ".txt";

	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_initial_gap_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_initial_gap_per_algo_cplex_cuts(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_gap_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_primal_gap_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_dual_gap_per_algo(algorithms.size(), std::vector<double>());
	std::vector<std::vector<double>> total_dual_gap_per_algo_no_cplex_cuts(algorithms.size(), std::vector<double>());

	std::vector<double> total_avg_initial_gap(algorithms.size(), 0.0);
	std::vector<double> total_avg_initial_gap_cplex_cuts(algorithms.size(), 0.0);
	std::vector<double> total_avg_gap(algorithms.size(), 0.0);
	std::vector<double> total_avg_primal_gap(algorithms.size(), 0.0);
	std::vector<double> total_avg_dual_gap(algorithms.size(), 0.0);
	std::vector<double> total_avg_dual_gap_no_cplex_cuts(algorithms.size(), 0.0);

	for (size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> initial_gap_per_algo(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> initial_gap_per_algo_cplex_cuts(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> gap_per_algo(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> primal_gap_per_algo(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> dual_gap_per_algo(algorithms.size(), std::vector<double>());
		std::vector<std::vector<double>> dual_gap_per_algo_no_cplex_cuts(algorithms.size(), std::vector<double>());

		std::vector<double> avg_initial_gap(algorithms.size(), 0.0);
		std::vector<double> avg_initial_gap_cplex_cuts(algorithms.size(), 0.0);
		std::vector<double> avg_gap(algorithms.size(), 0.0);
		std::vector<double> avg_primal_gap(algorithms.size(), 0.0);
		std::vector<double> avg_dual_gap(algorithms.size(), 0.0);
		std::vector<double> avg_dual_gap_no_cplex_cuts(algorithms.size(), 0.0);

		std::vector<double> st_dev_initial_gap(algorithms.size(), 0.0);
		std::vector<double> st_dev_initial_gap_cplex_cuts(algorithms.size(), 0.0);
		std::vector<double> st_dev_gap(algorithms.size(), 0.0);
		std::vector<double> st_dev_primal_gap(algorithms.size(), 0.0);
		std::vector<double> st_dev_dual_gap(algorithms.size(), 0.0);
		std::vector<double> st_dev_dual_gap_no_cplex_cuts(algorithms.size(), 0.0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);

		for (size_t i = 0; i < instances.size(); i++)
		{
			double root_lb = 0.0, root_ub2 = 0.0, root_ub = 0.0, best_lb = -1.0;
			double lb = 0.0, ub = 0.0;

			// find best primal bound
			for (size_t j = 0; j < algorithms.size(); j++)
			{
				std::string curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append((algorithms[j]).second);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lb1;
				std::string status;
				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				// std::cout << status << std::endl;
				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_lb1 << line.substr(pos + 2);

				if (s_lb1.str() == "-inf")
					lb = -1.0;
				else
					s_lb1 >> lb;

				if ((!double_equals(lb, -1.0)) && (status != "INFEASIBLE") && (double_greater(lb, best_lb)))
					best_lb = lb;

				// std::cout << "best_lb: " << best_lb << std::endl;
				input.close();
			}

			for (size_t j = 0; j < algorithms2.size(); j++)
			{
				double curr_initial_gap_cplex_cuts = 0.0, curr_initial_gap = 0.0, curr_gap = 0.0, curr_primal_gap = 0.0, curr_dual_gap = 0.0, curr_dual_gap_no_cplex_cuts = 0.0, time = 0.0;
				ub = 0.0;
				root_ub = 0.0;
				root_ub2 = 0.0;
				lb = 0.0;
				root_lb = 0.0;

				std::string curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append((algorithms2[j]).first);

				if (!user_cuts)
					curr_file += "_no_user_cuts";
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::string status, status2;
				std::string line;
				std::stringstream s_root_lb, s_root_ub2;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				getline(input, line);

				pos = line.find_first_of(":");
				s_root_lb << line.substr(pos + 2);
				if (s_root_lb.str() == "-inf")
					root_lb = -1.0;
				else
					s_root_lb >> root_lb;

				getline(input, line);

				pos = line.find_first_of(":");
				s_root_ub2 << line.substr(pos + 2);
				if (s_root_ub2.str() == "inf")
					root_ub2 = -1.0;
				else
					s_root_ub2 >> root_ub2;

				input.close();

				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append((algorithms2[j]).second);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_root_ub, s_lb, s_ub, s_time;
				getline(input, line);
				pos = line.find_first_of(":");
				status2 = line.substr(pos + 2);

				std::string::iterator end_pos2 = std::remove(status2.begin(), status2.end(), ' ');
				status2.erase(end_pos2, status2.end());

				getline(input, line);

				pos = line.find_first_of(":");
				s_root_ub << line.substr(pos + 2);
				if (s_root_ub.str() == "inf")
					root_ub = -1.0;
				else
					s_root_ub >> root_ub;

				getline(input, line);
				pos = line.find_first_of(":");
				s_lb << line.substr(pos + 2);
				if (s_lb.str() == "-inf")
					lb = -1.0;
				else
					s_lb >> lb;

				getline(input, line);
				pos = line.find_first_of(":");
				s_ub << line.substr(pos + 2);
				if (s_ub.str() == "inf")
					ub = -1.0;
				else
					s_ub >> ub;

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);
				s_time >> time;

				if ((status2 != "OPTIMAL") && (status2 != "INFEASIBLE"))
				{
					if ((double_equals(ub, -1)) || (double_equals(lb, -1)))
						curr_gap = 100.0;
					else if (double_equals(ub, 0.0) || double_less(time, time_limit))
						curr_gap = 0.0;
					else
						curr_gap = (100 * (ub - lb)) / ub;
				}
				else
				{
					curr_gap = 0.0;
				}

				if (status != "INFEASIBLE")
				{
					if ((double_equals(root_ub, -1)) || (double_equals(root_lb, -1)))
						curr_initial_gap = 100.0;
					else if (double_equals(root_ub, 0.0))
						curr_initial_gap = 0.0;
					else
						curr_initial_gap = (100 * (root_ub - root_lb)) / root_ub;

					if ((double_equals(root_ub2, -1)) || (double_equals(root_lb, -1)))
						curr_initial_gap_cplex_cuts = 100.0;
					else if (double_equals(root_ub2, 0.0))
						curr_initial_gap_cplex_cuts = 0.0;
					else
						curr_initial_gap_cplex_cuts = (100 * (root_ub2 - root_lb)) / root_ub2;
				}
				else
					curr_initial_gap = curr_initial_gap_cplex_cuts = 0.0;

				if (!double_equals(best_lb, -1.0))
				{
					if ((double_equals(root_lb, -1)) || (double_equals(root_ub2, -1)))
						curr_primal_gap = 100.0;
					else if (double_equals(root_ub2, 0.0))
						curr_primal_gap = 0.0;
					else
						curr_primal_gap = (100 * (best_lb - root_lb)) / root_ub2;

					if (double_equals(root_ub2, -1))
						curr_dual_gap = 100.0;
					else if (double_equals(root_ub2, 0.0))
						curr_dual_gap = 0.0;
					else
						curr_dual_gap = (100 * (root_ub2 - best_lb)) / root_ub2;

					if (double_equals(root_ub, -1))
						curr_dual_gap_no_cplex_cuts = 100.0;
					else if (double_equals(root_ub, 0.0))
						curr_dual_gap_no_cplex_cuts = 0.0;
					else
						curr_dual_gap_no_cplex_cuts = (100 * (root_ub - best_lb)) / root_ub;
				}
				else
					curr_primal_gap = curr_dual_gap = curr_dual_gap_no_cplex_cuts = 0.0;

				// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				initial_gap_per_algo[j].push_back(curr_initial_gap);
				initial_gap_per_algo_cplex_cuts[j].push_back(curr_initial_gap_cplex_cuts);
				gap_per_algo[j].push_back(curr_gap);
				primal_gap_per_algo[j].push_back(curr_primal_gap);
				dual_gap_per_algo[j].push_back(curr_dual_gap);
				dual_gap_per_algo_no_cplex_cuts[j].push_back(curr_dual_gap_no_cplex_cuts);

				total_initial_gap_per_algo[j].push_back(curr_initial_gap);
				total_initial_gap_per_algo_cplex_cuts[j].push_back(curr_initial_gap_cplex_cuts);
				total_gap_per_algo[j].push_back(curr_gap);
				total_primal_gap_per_algo[j].push_back(curr_primal_gap);
				total_dual_gap_per_algo[j].push_back(curr_dual_gap);
				total_dual_gap_per_algo_no_cplex_cuts[j].push_back(curr_dual_gap_no_cplex_cuts);

				total_avg_initial_gap[j] += curr_initial_gap;
				total_avg_initial_gap_cplex_cuts[j] += curr_initial_gap_cplex_cuts;
				total_avg_gap[j] += curr_gap;
				total_avg_primal_gap[j] += curr_primal_gap;
				total_avg_dual_gap[j] += curr_dual_gap;
				total_avg_dual_gap_no_cplex_cuts[j] += curr_dual_gap_no_cplex_cuts;

				avg_initial_gap[j] += curr_initial_gap;
				avg_initial_gap_cplex_cuts[j] += curr_initial_gap_cplex_cuts;
				avg_gap[j] += curr_gap;
				avg_primal_gap[j] += curr_primal_gap;
				avg_dual_gap[j] += curr_dual_gap;
				avg_dual_gap_no_cplex_cuts[j] += curr_dual_gap_no_cplex_cuts;

				input.close();
			}
			// getchar();getchar();
		}

		output << j + 1;
		if (stop)
			output << "\\_5\\%";

		for (size_t k = 0; k < algorithms2.size(); k++)
		{
			avg_initial_gap[k] /= (1.0 * (instances.size()));
			st_dev_initial_gap[k] = StDev(initial_gap_per_algo[k], avg_initial_gap[k]);
			avg_initial_gap_cplex_cuts[k] /= (1.0 * (instances.size()));
			st_dev_initial_gap_cplex_cuts[k] = StDev(initial_gap_per_algo_cplex_cuts[k], avg_initial_gap_cplex_cuts[k]);
			avg_gap[k] /= (1.0 * (instances.size()));
			st_dev_gap[k] = StDev(gap_per_algo[k], avg_gap[k]);
			avg_primal_gap[k] /= (1.0 * (instances.size()));
			st_dev_primal_gap[k] = StDev(primal_gap_per_algo[k], avg_primal_gap[k]);
			avg_dual_gap[k] /= (1.0 * (instances.size()));
			st_dev_dual_gap[k] = StDev(dual_gap_per_algo[k], avg_dual_gap[k]);
			avg_dual_gap_no_cplex_cuts[k] /= (1.0 * (instances.size()));
			st_dev_dual_gap_no_cplex_cuts[k] = StDev(dual_gap_per_algo_no_cplex_cuts[k], avg_dual_gap_no_cplex_cuts[k]);
			output << " & & " << avg_initial_gap_cplex_cuts[k] << " & " << st_dev_initial_gap_cplex_cuts[k] << " & & " << avg_primal_gap[k] << " & " << st_dev_primal_gap[k] << " & & " << avg_dual_gap[k] << " & " << st_dev_dual_gap[k];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for (size_t j = 0; j < algorithms2.size(); j++)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_initial_gap[j] /= (1.0 * ((total_initial_gap_per_algo[j]).size()));
		total_avg_initial_gap_cplex_cuts[j] /= (1.0 * ((total_initial_gap_per_algo_cplex_cuts[j]).size()));
		total_avg_gap[j] /= (1.0 * ((total_gap_per_algo[j]).size()));
		total_avg_primal_gap[j] /= (1.0 * ((total_primal_gap_per_algo[j]).size()));
		total_avg_dual_gap[j] /= (1.0 * ((total_dual_gap_per_algo[j]).size()));
		total_avg_dual_gap_no_cplex_cuts[j] /= (1.0 * ((total_dual_gap_per_algo_no_cplex_cuts[j]).size()));
		output << "& & " << total_avg_initial_gap_cplex_cuts[j] << " & " << StDev(total_initial_gap_per_algo_cplex_cuts[j], total_avg_initial_gap_cplex_cuts[j]) << "& & " << total_avg_primal_gap[j] << " & " << StDev(total_primal_gap_per_algo[j], total_avg_primal_gap[j]) << " & & " << total_avg_dual_gap[j] << " & " << StDev(total_dual_gap_per_algo[j], total_avg_dual_gap[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateLPImprovementsLatexTable(std::vector<std::string> dirs, bool stop)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_FBs_cb_csc3");
	algorithms.push_back("relax_FBs_GCCs_cb_csc3");
	algorithms.push_back("relax_FBs_CCs_cb_csc3");
	// algorithms.push_back("relax_FBs_CCCs_cb_csc3");
	algorithms.push_back("relax_FBs_LCIs_cb_csc3");
	// algorithms.push_back("relax_FBs_AVICs_cb_csc3");

	algorithms.push_back("relax_FBs_GCCs_CCs_cb_csc3");
	algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_cb_csc3");
	// algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_AVICs_cb_csc3");
	// algorithms.push_back("relax_FBs_CCCs_LCIs_cb_csc3");
	// algorithms.push_back("relax_FBs_CCCs_LCIs_AVICs_cb_csc3");*/
	// algorithms.push_back("relax_5_cb_csc3");
	/*algorithms.push_back("relax_6_cb_csc3");
	algorithms.push_back("relax_7_cb_csc3");
	algorithms.push_back("relax_8_cb_csc3");
	algorithms.push_back("relax_9_cb_csc3");*/

	// algorithms.push_back("cb_csc3");
	// algorithms.push_back("cb_angle_0.03_clique_ccs_lcis_csc3");
	/*algorithms.push_back("relax_all_gccs");
	//algorithms.push_back("relax_all_gccs2");
	//algorithms.push_back("relax_all_gccs3");
	//algorithms.push_back("relax_gccs_cb_csc3");
	algorithms.push_back("relax_all_ccs");
	//algorithms.push_back("relax_ccs_cb_csc3");
	algorithms.push_back("relax_all_gccs_ccs");
	//algorithms.push_back("relax_gccs_ccs_cb_csc3");
	//algorithms.push_back("relax_all_gccs_ccs_lcis");
	//algorithms.push_back("relax_clique_ccs");
	//algorithms.push_back("relax_clique_ccs_active_vertices");
	algorithms.push_back("relax_clique_ccs_active_vertices_unilateral");
	//algorithms.push_back("relax_clique_ccs_total_active_vertices_unilateral");
	algorithms.push_back("relax_maximum_clique_ccs_active_vertices_unilateral");
	algorithms.push_back("relax_angle_0.1_maximum_clique_ccs_active_vertices_unilateral");
	algorithms.push_back("relax_angle_0.05_clique_ccs_active_vertices_unilateral");*/
	// algorithms.push_back("relax_csc3e");
	// algorithms.push_back("relax_cb_csc3e"); // é o relax_clique_ccs_cb_csc3e
	// algorithms.push_back("relax_clique_ccs_lcis_cb_csc3e");
	// algorithms.push_back("relax_angle_0.03_clique_ccs_active_vertices_unilateral");
	// algorithms.push_back("relax_lcis_clique_ccs_active_vertices_unilateral");
	// algorithms.push_back("relax_clique_ccs");
	// algorithms.push_back("relax_all_cb_csc3");
	// algorithms.push_back("relax_all_maximal_cliques_ccs_cb_csc3");
	// algorithms.push_back("relax_all_clique_ccs_cb_csc3");
	// algorithms.push_back("relax_all_clique_ccs_cb_csc3_new");
	// algorithms.push_back("relax_clique_ccs_cb_csc3");
	// algorithms.push_back("relax_clique_ccs_lcis_cb_csc3");
	// algorithms.push_back("relax_clique_ccs_cb_csc3_new_unilateral");
	// algorithms.push_back("relax_clique_ccs_cb_csc3");
	/*algorithms.push_back("relax_gccs_cb_csc3");
	algorithms.push_back("relax_ccs_cb_csc3");
	algorithms.push_back("relax_lcis_cb_csc3");
	algorithms.push_back("relax_gccs_ccs_cb_csc3");
	algorithms.push_back("relax_all_cb_csc3");*/
	// algorithms.push_back("relax_cccs_cb_csc3");
	// algorithms.push_back("relax_all_maximal_cliques_ccs_cb_csc3");
	// algorithms.push_back("relax_gccs_all_maximal_cliques_ccs_cb_csc3");

	// algorithms.push_back("relax_3gccs_cb_csc3");

	std::fstream output;
	std::string output_name;
	if (stop)
		output_name = ".//tables//latex//stop_table_LP_improvements.txt";
	else
		output_name = ".//tables//latex//table_LP_improvements.txt";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;

	std::vector<std::vector<double>> total_improvement_per_algo(algorithms.size(), std::vector<double>());
	std::vector<double> total_avg_improvement(algorithms.size(), 0.0);

	for (size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::vector<double>> improvement_per_algo(algorithms.size(), std::vector<double>());
		std::vector<double> avg_improvement(algorithms.size(), 0.0);
		std::vector<double> st_dev(algorithms.size(), 0.0);

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);
		for (size_t i = 0; i < instances.size(); i++)
		{
			std::cout << instances[i] << std::endl;
			double original_lp = 0.0;
			for (size_t j = 0; j < algorithms.size(); j++)
			{
				double curr_improvement = 0.0;
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp;
				std::string status;
				double lp = 0.0;
				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);

				if (s_lp.str() == "inf")
					lp = -1;
				else
					s_lp >> lp;

				if (j == 0)
				{
					if ((double_equals(lp, -1)) || (status == "INFEASIBLE"))
					{
						original_lp = -1.0;
					}
					else
						original_lp = lp;
					curr_improvement = 0.0;
				}
				else
				{
					if ((double_equals(lp, -1)) || (status == "INFEASIBLE") || (double_equals(original_lp, -1)))
					{
						curr_improvement = -1;
					}
					else
					{
						if (double_equals(original_lp, 0.0))
							curr_improvement = 0.0;
						else
							curr_improvement = (100 * (original_lp - lp)) / original_lp;

						if (double_less(curr_improvement, 0))
							std::cout << original_lp << " - " << lp << std::endl;
					}
				}

				// std::cout << "lp: " << lp << " improvement: " << curr_improvement << std::endl;

				if (!double_equals(curr_improvement, -1))
				{
					improvement_per_algo[j].push_back(curr_improvement);
					total_improvement_per_algo[j].push_back(curr_improvement);
					total_avg_improvement[j] += curr_improvement;
					avg_improvement[j] += curr_improvement;
				}
				input.close();
			}

			// getchar();getchar();
		}

		output << j + 1;

		for (size_t j = 1; j < algorithms.size(); j++)
		{
			avg_improvement[j] /= (1.0 * ((improvement_per_algo[j]).size()));
			st_dev[j] = StDev(improvement_per_algo[j], avg_improvement[j]);
			output << " & & " << avg_improvement[j] << " & " << st_dev[j];
		}
		output << "\\\\" << std::endl;
	}

	output << "Total";
	for (size_t j = 1; j < algorithms.size(); j++)
	{
		// std::cout << total_improvement_per_algo[j].size() << std::endl;
		total_avg_improvement[j] /= (1.0 * ((total_improvement_per_algo[j]).size()));
		output << "& & " << total_avg_improvement[j] << " & " << StDev(total_improvement_per_algo[j], total_avg_improvement[j]);
	}

	output << "\\\\" << std::endl;

	output.close();
}

void GenerateCutsConfigurationsPerInstanceResultsCSV(std::vector<std::string> dirs, bool stop)
{
	std::string curr_file;
	std::vector<std::string> algorithms;

	algorithms.push_back("relax_FBs_cb_csc3");
	algorithms.push_back("relax_FBs_GCCs_cb_csc3");
	algorithms.push_back("relax_FBs_CCs_cb_csc3");
	// algorithms.push_back("relax_FBs_CCCs_cb_csc3");
	algorithms.push_back("relax_FBs_LCIs_cb_csc3");
	// algorithms.push_back("relax_FBs_AVICs_cb_csc3");

	algorithms.push_back("relax_FBs_GCCs_CCs_cb_csc3");
	algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_cb_csc3");
	// algorithms.push_back("relax_FBs_GCCs_CCs_LCIs_AVICs_cb_csc3");
	// algorithms.push_back("relax_FBs_CCCs_LCIs_cb_csc3");
	// algorithms.push_back("relax_FBs_CCCs_LCIs_AVICs_cb_csc3");

	std::fstream output;
	std::string output_name;
	if (stop)
		output_name = ".//tables//CSV//stop_table_LP_bounds.csv";
	else
		output_name = ".//tables//CSV//table_LP_bounds.csv";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << "Configuration;";

	for (size_t i = 0; i < algorithms.size(); ++i)
	{
		output << i << ";";
	}

	output << std::endl;

	output << "Instance;";

	for (size_t i = 0; i < algorithms.size(); ++i)
		output << "upper bound;";

	output << std::endl;
	output << std::setprecision(2) << std::fixed;

	for (size_t j = 0; j < dirs.size(); j++)
	{

		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);

		std::sort(instances.begin(), instances.end());

		for (size_t i = 0; i < instances.size(); i++)
		{
			std::cout << instances[i] << std::endl;

			if (stop)
				output << instances[i].substr(0, instances[i].size() - 4) << "_5%;";
			else
				output << instances[i].substr(0, instances[i].size() - 4) << ";";

			for (size_t j = 0; j < algorithms.size(); j++)
			{
				double curr_improvement = 0.0;
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[j]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp;
				std::string status;
				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);

				if (status == "INFEASIBLE")
					output << " - ;";
				else
					output << s_lp.str() << ";";

				// output << s_lp.str() << ";";
				input.close();
			}

			// getchar();getchar();
			output << std::endl;
		}
	}

	output.close();
}

void GenerateExactAlgorithmsPerInstanceResultsCSV(std::vector<std::string> dirs, bool stop, double time_limit)
{
	std::string curr_file;
	double time = 0.0;
	double ub = 0.0, lb = 0.0;

	std::vector<std::string> algorithms;
	algorithms.push_back("bc_b1");
	algorithms.push_back("cb_csc3");
	/*algorithms.push_back("cb2_csc3");*/

	/*std::string algo = "cb2_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	algo = "cb2_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);*/

	std::fstream output;
	std::string output_name;

	// std::cout << output_name << std::endl;

	if (stop)
		output_name = ".//tables//CSV//stop_exact_algorithms_per_instance.csv";
	else
		output_name = ".//tables//CSV//exact_algorithms_per_instance.csv";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	output << std::setprecision(2) << std::fixed;

	output << "Algorithm;";
	for (size_t i = 0; i < algorithms.size(); ++i)
	{
		for (size_t j = 0; j < 3; ++j)
		{
			output << algorithms[i] << ";";
		}
	}

	output << std::endl;

	output << "Instance;";
	for (size_t i = 0; i < algorithms.size(); ++i)
	{
		output << "lb;ub;time(s);";
	}

	output << std::endl;

	for (size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);

		std::sort(instances.begin(), instances.end());

		// output << folder << std::endl;

		for (size_t i = 0; i < instances.size(); i++)
		{
			if (stop)
				output << instances[i].substr(0, instances[i].size() - 4) << "_5%;";
			else
				output << instances[i].substr(0, instances[i].size() - 4) << ";";

			for (size_t k = 0; k < algorithms.size(); ++k)
			{
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[k]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lb, s_ub, s_time;
				std::string status;
				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_lb << line.substr(pos + 2);
				getline(input, line);
				pos = line.find_first_of(":");
				s_ub << line.substr(pos + 2);

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);

				s_time >> time;
				s_lb >> lb;
				s_ub >> ub;
				// gap1 = 0.0;

				// if((status != "OPTIMAL") && (status != "INFEASIBLE"))
				//{
				// if(!double_less(time,time_limit)) gap1 = (100.0*(ub-lb))/ub;
				//}

				if (((status == "OPTIMAL") || (double_less(time, 7200))) && (!double_equals(lb, ub)))
					std::cout << algorithms[k] << "_" << instances[i] << " " << status << " " << time << "s " << lb << " != " << ub << std::endl;

				if (double_greater(time, 7200))
					time = 7200;
				if (status == "INFEASIBLE")
					output << " - ; - ; " << time << ";";
				else
					output << s_lb.str() << " ; " << s_ub.str() << " ; " << time << ";";
				input.close();
			}
			output << std::endl;
		}
	}
	output.close();
}

void GenerateAppendixLatexTable(std::vector<std::string> dirs, bool stop, double time_limit)
{
	std::string curr_file;
	double time = 0.0;
	double ub = 0.0, lb = 0.0;

	std::vector<std::string> algorithms;
	/*algorithms.push_back("bc_b1");
	algorithms.push_back("cb_csc3");
	algorithms.push_back("cb2_csc3");*/

	std::string algo = "cb2_csc3_";
	if (stop)
		algo += "stop_";
	algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	algo = "cb2_csc3_";
	if (stop)
		algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	std::fstream output;
	std::string output_name;

	// std::cout << output_name << std::endl;

	if (stop)
		output_name = ".//tables//latex//stop_appendix.txt";
	else
		output_name = ".//tables//latex//appendix.txt";
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	output << std::setprecision(2) << std::fixed;

	for (size_t j = 0; j < dirs.size(); j++)
	{
		std::vector<std::string> instances;
		std::string folder = dirs[j].substr(20);
		AddInstancesFromDirectory(dirs[j], instances, false);

		std::sort(instances.begin(), instances.end());

		// output << folder << std::endl;

		for (size_t i = 0; i < instances.size(); i++)
		{
			if (stop)
				output << instances[i].substr(0, instances[i].size() - 4) << "\\_5\\%";
			else
				output << instances[i].substr(0, instances[i].size() - 4);

			for (size_t k = 0; k < algorithms.size(); ++k)
			{
				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if (stop)
					curr_file.append("stop_");
				curr_file.append(algorithms[k]);
				curr_file.append("_");
				curr_file.append(instances[i]);
				// std::cout << curr_file << std::endl;

				std::fstream input;
				input.open(curr_file.c_str(), std::fstream::in);

				if (!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lb, s_ub, s_time;
				std::string status;
				std::string line;

				getline(input, line);
				size_t pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
				status.erase(end_pos, status.end());

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_lb << line.substr(pos + 2);
				getline(input, line);
				pos = line.find_first_of(":");
				s_ub << line.substr(pos + 2);

				getline(input, line);
				getline(input, line);
				pos = line.find_first_of(":");
				s_time << line.substr(pos + 2);

				s_time >> time;
				s_lb >> lb;
				s_ub >> ub;
				// gap1 = 0.0;

				// if((status != "OPTIMAL") && (status != "INFEASIBLE"))
				//{
				// if(!double_less(time,time_limit)) gap1 = (100.0*(ub-lb))/ub;
				//}

				if (((status == "OPTIMAL") || (double_less(time, 7200))) && (!double_equals(lb, ub)))
					std::cout << algorithms[k] << "_" << instances[i] << " " << status << " " << time << "s " << lb << " != " << ub << std::endl;

				if (k != 0)
					output << " & ";
				if (double_greater(time, 7200))
					time = 7200;
				if (status == "INFEASIBLE")
					output << " & -- & -- & " << time;
				else
					output << " & " << s_lb.str() << " & " << s_ub.str() << " & " << time;
				input.close();

				/*curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append("relax_none_cb_csc3_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				input.open(curr_file.c_str(),std::fstream::in);
				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lp, s_lp_cuts;

				getline(input,line);
				getline(input,line);
				pos = line.find_first_of(":");
				s_lp << line.substr(pos + 2);

				input.close();

				curr_file = ".//solutions//";
				curr_file.append(folder);
				curr_file.append("s_");
				if(stop) curr_file.append("stop_");
				curr_file.append("cb_csc3_");
				curr_file.append(instances[i]);
				//std::cout << curr_file << std::endl;

				input.open(curr_file.c_str(),std::fstream::in);
				if(!input.is_open())
				{
					std::cout << "Could not open file " << curr_file << std::endl;
					throw 4;
				}

				std::stringstream s_lb2, s_ub2, s_time2;

				getline(input,line);
				pos = line.find_first_of(":");
				status = line.substr(pos + 2);

				end_pos = std::remove(status.begin(),status.end(),' ');
				status.erase(end_pos, status.end());

				getline(input,line);
				pos = line.find_first_of(":");
				s_lp_cuts << line.substr(pos + 2);
				getline(input,line);
				pos = line.find_first_of(":");
				s_lb2 << line.substr(pos + 2);
				getline(input,line);
				pos = line.find_first_of(":");
				s_ub2 << line.substr(pos + 2);

				getline(input,line);
				getline(input,line);
				pos = line.find_first_of(":");
				s_time2 << line.substr(pos + 2);

				s_time2 >> time;

				s_lb2 >> lb;
				s_ub2 >> ub;
				gap2 = 0.0;

				if((status != "OPTIMAL") && (status != "INFEASIBLE"))
				{
					if(!double_less(time,time_limit)) gap2 = (100.0*(ub-lb))/ub;
				}

				if(double_greater(time,7200)) time = 7200;
				if(status == "INFEASIBLE") output << "-- & -- & " << time << " & -- & --";
				else output << s_lb2.str() << " & " << s_ub2.str() << " & " << time << " & " << s_lp.str() << " & " << s_lp_cuts.str();*/
				// for(int i = 0; i < 10; i++) getline(input,line);
			}
			output << "\\\\" << std::endl;
		}

		if (j != dirs.size() - 1)
			output << "\\\\" << std::endl;
	}
	output.close();
}

void GenerateCSVFile(std::string dir, double time_limit, bool stop)
{
	std::string folder = dir.substr(20);

	std::vector<std::string> instances;
	AddInstancesFromDirectory(dir, instances, false);
	std::sort(instances.begin(), instances.end());
	std::string curr_file;

	std::vector<std::string> algorithms;
	/*algorithms.push_back("bcsc1");
	algorithms.push_back("bcsc2");
	algorithms.push_back("bcsc3");*/
	// algorithms.push_back("bcsc4");

	/*algorithms.push_back("bcmc1");
	algorithms.push_back("bcmc2");
	algorithms.push_back("bcmc3");*/
	// algorithms.push_back("bcmc4");

	// algorithms.push_back("bctc");
	// algorithms.push_back("bcmtc");

	// algorithms.push_back("csc1");
	// algorithms.push_back("csc3");
	// algorithms.push_back("csc4");

	// algorithms.push_back("cmc1");
	// algorithms.push_back("cmc2");

	// algorithms.push_back("ctc1");
	// algorithms.push_back("ctc2");

	// algorithms.push_back("cmtc1");
	// algorithms.push_back("cmtc2");

	// algorithms.push_back("kb");

	// algorithms.push_back("b1");
	// algorithms.push_back("b2");

	/*algorithms.push_back("bcsc4");
	algorithms.push_back("bcmc4");
	algorithms.push_back("bctc");
	algorithms.push_back("bcmtc");*/
	// algorithms.push_back("csc3_relax1");

	/*algorithms.push_back("csc3_raw");
	algorithms.push_back("csc3");
	algorithms.push_back("csc3_avics");
	algorithms.push_back("avics_csc3");*/

	// algorithms.push_back("bc_csc3_heuristic_all_cuts");
	// algorithms.push_back("csc3");
	// algorithms.push_back("csc3_heuristic_advance_2");
	// algorithms.push_back("cb_csc3_heuristic_avics_hard_coded_advance_2");
	// algorithms.push_back("cb_csc3_heuristic_avics_cuts_advance_2");

	// algorithms.push_back("cb_csc3");
	// algorithms.push_back("cb_csc3_heuristic");
	// algorithms.push_back("cb_csc3_heuristic_advance_2");
	// algorithms.push_back("cb_csc3_new_heuristic_advance_2");

	/*algorithms.push_back("cb_avics_gccs_ccs_lcis_csc3");
	algorithms.push_back("cb_avics_clique_ccs_lcis_csc3");
	algorithms.push_back("cb_angle_0.03_clique_ccs_lcis_csc3");*/
	// algorithms.push_back("cb_avics_clique_ccs_lcis_csc3");
	// algorithms.push_back("cb_heuristic");
	// algorithms.push_back("cb_heuristic_advance_2");
	// algorithms.push_back("cb_csc3_heuristic_avics_cover_advance_2");
	algorithms.push_back("bc_b1");
	algorithms.push_back("cb_csc3");
	// algorithms.push_back("cb2_csc3");
	// algorithms.push_back("bc_b1");

	std::string algo = "cb_csc3_";
	/*if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	algo = "cb_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	/*algo = "cb2_csc3_";
	if(stop) algo += "stop_";
	algo += "fp_cuts_alpha2_0_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);*/

	algo = "cb2_csc3_";
	if (stop)
		algo += "stop_";
	algo += "fp_cuts_alpha2_100_p_fixed_10_PR_ALNS_5000";
	algorithms.push_back(algo);

	// algorithms.push_back("cb_angle_0.03_clique_ccs_lcis_csc3");
	/*algorithms.push_back("bc_branchcallback_csc3");
	algorithms.push_back("bc_cutcallback_branchcallback_csc3");
	algorithms.push_back("bc_cutcallback_csc3");
	algorithms.push_back("bc_initialflowbounds_cutcallback_branchcallback_csc3");*/
	// algorithms.push_back("st_bc_b1");
	// algorithms.push_back("st_cb_csc3");
	// algorithms.push_back("csc3_relax3");
	// algorithms.push_back("cb_r3");
	// algorithms.push_back("cb_r4");
	/*algorithms.push_back("bc1_0.05_PREPROCESS_WITH_Y_b1");
	algorithms.push_back("bc1_0.05_PREPROCESS_WITH_Y_csc3");
	algorithms.push_back("bc4_0.6_0.1_PREPROCESS_WITH_Y_b1");
	algorithms.push_back("bc4_0.6_0.1_PREPROCESS_WITH_Y_csc3");*/
	// algorithms.push_back("csc4");
	/*algorithms.push_back("cmc2");
	algorithms.push_back("ctc1");
	algorithms.push_back("cmtc1");
	algorithms.push_back("kb");*/
	// algorithms.push_back("b2");

	std::fstream output;
	std::string output_name = ".//tables//CSV//";
	output_name.append(folder);
	if (stop)
		output_name.append("correct_table_stop.csv");
	else
		output_name.append("correct_table_top.csv");
	output.open(output_name.c_str(), std::fstream::out);

	if (!output.is_open())
	{
		std::cout << "Could not open file " << output_name << std::endl;
		throw 1;
	}

	// std::cout << output_name << std::endl;

	output << std::setprecision(2) << std::fixed;
	output << "Formulation;;;;";

	for (size_t j = 0; j < algorithms.size(); j++)
	{
		if (stop)
			algorithms[j] = "stop_" + algorithms[j];
		output << algorithms[j] << ";" << algorithms[j] << ";" << algorithms[j] << ";" << algorithms[j] << ";" << algorithms[j] << ";" << algorithms[j] << ";;";
	}

	output << std::endl
		   << "Instance;T;m;;";
	for (size_t j = 0; j < algorithms.size(); j++)
	{
		output << "status;lp;lb;ub;gap;time(s);;";
	}

	output << std::endl;

	for (size_t i = 0; i < instances.size(); i++)
	{
		curr_file = dir + instances[i];
		double limit = 0.0;
		int num_vehicles = 0;

		std::fstream file;
		file.open(curr_file.c_str(), std::fstream::in);

		if (file.is_open())
		{
			std::stringstream parameter;
			std::string line;
			getline(file, line);

			getline(file, line);

			getline(file, line);
			std::size_t pos_1 = line.find_first_of(" ");

			parameter.clear();
			parameter << line.substr(pos_1 + 1);
			parameter >> num_vehicles;

			getline(file, line);
			pos_1 = line.find_first_of(" ");

			parameter.clear();
			parameter << line.substr(pos_1 + 1);
			parameter >> limit;

			file.close();
		}

		output << instances[i] << ";" << limit << ";" << num_vehicles << ";;";
		for (size_t j = 0; j < algorithms.size(); j++)
		{
			curr_file = ".//solutions//";
			curr_file.append(folder);
			curr_file.append("s_");
			curr_file.append(algorithms[j]);
			curr_file.append("_");
			curr_file.append(instances[i]);
			// std::cout << curr_file << std::endl;

			std::fstream input;
			std::cout << curr_file << std::endl;
			input.open(curr_file.c_str(), std::fstream::in);

			if (!input.is_open())
			{
				output << ";"
					   << ";"
					   << ";"
					   << ";"
					   << ";"
					   << ";;";
				continue;
			}

			std::stringstream s_lp, s_lb, s_ub, s_time;
			std::string status;
			double lp, lb, ub, time, gap;
			std::string line;

			getline(input, line);
			size_t pos = line.find_first_of(":");
			status = line.substr(pos + 2);

			getline(input, line);
			pos = line.find_first_of(":");
			s_lp << line.substr(pos + 2);

			getline(input, line);
			pos = line.find_first_of(":");
			s_lb << line.substr(pos + 2);

			getline(input, line);
			pos = line.find_first_of(":");
			s_ub << line.substr(pos + 2);

			getline(input, line);
			getline(input, line);
			pos = line.find_first_of(":");
			s_time << line.substr(pos + 2);

			s_lp >> lp;
			s_lb >> lb;
			s_ub >> ub;
			s_time >> time;

			std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
			status.erase(end_pos, status.end());

			if (double_equals(lp, -1))
				output << status << ";"
					   << "-"
					   << ";";
			else
				output << status << ";" << lp << ";";
			if (double_less(time, time_limit))
				status = "OPTIMAL";

			if (status != "INFEASIBLE")
			{
				gap = (100.0 * (ub - lb)) / ub;

				output << lb << ";" << ub << ";" << gap << ";" << time << ";;";
			}
			else
				output << "-"
					   << ";"
					   << "-"
					   << ";"
					   << "-"
					   << ";" << time << ";;";
		}
		output << std::endl;
	}

	output.close();
}

double RetrieveTimeSpentInFP(Instance &inst, std::string algo, std::string folder, std::string file_name)
{
	double t1 = 0.0, t2 = 0.0;
	std::stringstream s_t1, s_t2;

	std::fstream input;
	std::string path = ".//solutions";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	input.open(path.c_str(), std::fstream::in);

	if (!(input.is_open()))
		throw "Error opening file 2";
	std::string line;
	for (int i = 1; i <= 7; ++i)
		getline(input, line);

	size_t pos = line.find_first_of(":");
	s_t1 << line.substr(pos + 2);
	s_t1 >> t1;

	for (int i = 1; i <= 5; ++i)
		getline(input, line);

	pos = line.find_first_of(":");
	s_t2 << line.substr(pos + 2);
	s_t2 >> t2;

	input.close();

	// std::cout << path << std::endl;
	// std::cout << t1 << std::endl;
	// std::cout << t2 << std::endl;
	return t1 + t2;
}

double RetrieveTimeSpentInLNS(Instance &inst, std::string algo, std::string folder, std::string file_name)
{
	double time = 0.0;
	std::stringstream s_time;

	std::fstream input;
	std::string path = ".//solutions";
	path.append(folder);
	// struct stat sb;
	// if(stat(path.c_str(),&sb) != 0 || !S_ISDIR(sb.st_mode)) mkdir(path.c_str(),0777);
	path.append("s_");
	path.append(algo);
	path.append("_");
	path.append(file_name);
	// std::cout << path << std::endl;

	input.open(path.c_str(), std::fstream::in);

	if (!(input.is_open()))
		throw "Error opening file 2";
	std::string line;
	for (int i = 1; i <= 5; ++i)
		getline(input, line);

	size_t pos = line.find_first_of(":");
	s_time << line.substr(pos + 2);
	s_time >> time;
	input.close();

	// std::cout << path << std::endl;
	// std::cout << time << std::endl;
	return time;
}

void AddInstances(std::vector<std::string> &instances)
{
	// AddInstancesFromDirectory(".//instances//STOP//unsolved//",instances);
	AddInstancesFromDirectory(".//instances//STOP//current//", instances);
	/*AddInstancesFromDirectory(".//instances//STOP//Set_21_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_32_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_33_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_64_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_66_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_100_234//Set_0.05//",instances);
	AddInstancesFromDirectory(".//instances//STOP//Set_102_234//Set_0.05//",instances);*/
}

static const struct option longOpts[] = {
	{"solution-dir", required_argument, NULL, 'a'},
	{"compact", no_argument, NULL, 'b'},
	{"cutting-plane", no_argument, NULL, 'c'},
	{"branch-and-cut", no_argument, NULL, 'd'},
	{"baseline", no_argument, NULL, 'e'},
	{"capacity-based", no_argument, NULL, 'f'},
	{"time-limit", required_argument, NULL, 'g'},
	{"instance", required_argument, NULL, 'h'},
	{"generate-convex-hull", no_argument, NULL, 'i'},
	{"num-vehicles", required_argument, NULL, 'j'},
	{"service-time-deviation", required_argument, NULL, 'k'},
	{"uncertainty-budget", required_argument, NULL, 'l'},
	{"CCCs", no_argument, NULL, 'm'},
	{"AVICs", no_argument, NULL, 'n'},
	{"solve-relaxed", no_argument, NULL, 'o'},
	{"solve-benders", no_argument, NULL, 'p'},
	{"benders-generic-callback", no_argument, NULL, 'q'},
	{"combine-feas-opt-cuts", no_argument, NULL, 'r'},
	{"separate-benders-cuts-relaxation", no_argument, NULL, 's'},
	{"warm-start", no_argument, NULL, 't'},
	{"kernel-search", no_argument, NULL, 'u'},
	{"ks-max-size-bucket", required_argument, NULL, 'v'},
	{"ks-min-time-limit", required_argument, NULL, 'w'},
	{"ks-max-time-limit", required_argument, NULL, 'x'},
	{"ks-decay-factor", required_argument, NULL, 'y'},
	{"ks-feasibility-emphasis", required_argument, NULL, 'z'},
	{"seed", required_argument, NULL, 'A'},
	{"alns", no_argument, NULL, 'B'},
	{"alns-num-iterations", required_argument, NULL, 'C'},
	{"alns-pool-size", required_argument, NULL, 'D'},
	{"metaheuristic-multi-threading", no_argument, NULL, 'E'}, // Exact and Kernel search will always be multithreading!
	{"simulated-annealing", no_argument, NULL, 'F'},
	{"sa-temperature-decay-rate", required_argument, NULL, 'G'},
	{"sa-sampling-size", required_argument, NULL, 'H'},
	{"sa-target-acceptance-probability", required_argument, NULL, 'I'},
	{NULL, no_argument, NULL, 0}};

void ParseArgumentsAndRun(int argc, char *argv[])
{
	std::string instance, folder, file_name, dir_solutions;
	int c;
	int seed = 0;
	int num_vehicles = 0, uncertainty_budget = 0;
	bool solve_compact = false, solve_cb = false, solve_bc = false, baseline = false, capacity_based = false, generate_convex_hull = false;
	bool solve_relaxed = false, solve_benders = false, solve_generic_callback = false, combine_feas_opt_cuts = false;
	double time_limit = -1.0, service_time_deviation = 0.0;
	bool force_use_all_vehicles = false;
	bool export_model = false;
	bool separate_benders_cuts_relaxation = false;
	bool compute_initial_solution_heuristic = false;
	bool solve_kernel_search = false;
	bool solve_alns = false;
	bool metaheuristic_multi_threading = false;
	int ks_max_size_bucket = K_KS_MAX_SIZE_BUCKET, ks_min_time_limit = K_KS_MIN_TIME_LIMIT, ks_max_time_limit = K_KS_MAX_TIME_LIMIT;
	double ks_decay_factor = K_KS_DECAY_FACTOR_TIME_LIMIT;
	bool ks_feasibility_emphasis = false;
	bool solve_simulated_annealing = false;
	HeuristicSolution *initial_sol = nullptr;
	int alns_num_iterations = K_ALNS_NUM_ITERATIONS, alns_pool_size = K_ALNS_SIZE_OF_POOL;
	double sa_temperature_decay_rate = K_SA_TEMP_DECAY_RATE, sa_target_acceptance_probability = K_SA_TARGET_ACCEPTANCE_PROBABILITY;
	size_t sa_sampling_size = K_SA_SAMPLING_SIZE;

	std::vector<bool> *CALLBACKS_SELECTION = GetCallbackSelection();

	while ((c = getopt_long(argc, argv, "g:h:j:k:l:v:w:x:y:z:A:C:D:G:H:I:", longOpts, NULL)) != -1)
	{
		switch (c)
		{
		case 'a':
			dir_solutions = std::string(optarg);
			break;
		case 'b':
			solve_compact = true;
			break;
		case 'c':
			solve_cb = true;
			break;
		case 'd':
			solve_bc = true;
			break;
		case 'e':
			baseline = true;
			break;
		case 'f':
			capacity_based = true;
			break;
		case 'g':
			if (optarg)
				time_limit = std::atoi(optarg);
			break;
		case 'h':
			if (optarg)
				instance = std::string(optarg);
			break;
		case 'i':
			generate_convex_hull = true;
			break;
		case 'j':
			if (optarg)
				num_vehicles = std::atoi(optarg);
			break;
		case 'k':
			if (optarg)
				service_time_deviation = std::atof(optarg);
			break;
		case 'l':
			if (optarg)
				uncertainty_budget = std::atoi(optarg);
			break;
		case 'm':
			(*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT] = true;
			break;
		case 'n':
			(*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT] = true;
			break;
		case 'o':
			solve_relaxed = true;
			break;
		case 'p':
			solve_benders = true;
			break;
		case 'q':
			solve_generic_callback = true;
			break;
		case 'r':
			combine_feas_opt_cuts = true;
			break;
		case 's':
			separate_benders_cuts_relaxation = true;
			break;
		case 't':
			compute_initial_solution_heuristic = true;
			break;
		case 'u':
			solve_kernel_search = true;
			break;
		case 'v':
			ks_max_size_bucket = std::atoi(optarg);
			break;
		case 'w':
			ks_min_time_limit = std::atoi(optarg);
			break;
		case 'x':
			ks_max_time_limit = std::atoi(optarg);
			break;
		case 'y':
			ks_decay_factor = std::atof(optarg);
			break;
		case 'z':
			ks_feasibility_emphasis = (std::atoi(optarg) != 0);
			break;
		case 'A':
			if (optarg)
				seed = std::atoi(optarg);
			break;
		case 'B':
			solve_alns = true;
			break;
		case 'C':
			if (optarg)
				alns_num_iterations = std::atoi(optarg);
			break;
		case 'D':
			if (optarg)
				alns_pool_size = std::atoi(optarg);
			break;
		case 'E':
			metaheuristic_multi_threading = true;
			break;
		case 'F':
			solve_simulated_annealing = true;
			break;
		case 'G':
			if (optarg)
				sa_temperature_decay_rate = std::atof(optarg);
			break;
		case 'H':
			if (optarg)
				sa_sampling_size = std::atoi(optarg);
			break;
		case 'I':
			if (optarg)
				sa_target_acceptance_probability = std::atof(optarg);
			break;
		}
	}

	srand(seed);

	if (((solve_compact || solve_cb || solve_benders) && !baseline && !capacity_based) ||
		(solve_compact && solve_bc) || (solve_compact && solve_benders) ||
		(solve_compact && solve_cb) || (solve_bc && solve_cb) ||
		(solve_compact && solve_kernel_search) || (solve_cb && solve_kernel_search) ||
		(solve_relaxed && compute_initial_solution_heuristic) || (solve_relaxed && solve_alns) ||
		(solve_relaxed && solve_simulated_annealing) || (solve_relaxed && solve_kernel_search) ||
		(solve_benders && solve_kernel_search) || (solve_alns && solve_simulated_annealing))
		throw 2;
	if (baseline && capacity_based)
		throw 3;

	if (!dir_solutions.empty())
		dir_solutions += "//";

	split_file_path(instance, folder, file_name);
	// std::cout << "* " << folder << " " << file_name << std::endl;

	// std::cout << num_vehicles << " " << service_time_deviation << " " << uncertainty_budget << std::endl;
	Instance inst(folder, file_name, num_vehicles, service_time_deviation, uncertainty_budget, false);
	const Graph *graph = inst.graph();
	std::string instance_name = inst.GetInstanceName();
	std::cout << instance_name << std::endl;

	Solution<double> sol(graph->num_vertices());

	double *R0 = Dijkstra(graph, false, false);
	double *Rn = Dijkstra(graph, true, false);

	// std::cout << inst << std::endl;
	if (compute_initial_solution_heuristic)
	{
		if (solve_kernel_search)
		{
			Formulation formulation;
			if (baseline)
			{
				formulation = Formulation::baseline;
			}
			if (capacity_based)
			{
				formulation = Formulation::single_commodity;
			}

			KernelSearch ks(inst);
			std::cout << KSHeuristicSolution::GenerateFileName(formulation, ks_max_size_bucket, ks_min_time_limit, ks_max_time_limit, ks_decay_factor, ks_feasibility_emphasis) << std::endl;
			initial_sol = ks.Run(formulation, ks_max_size_bucket, ks_min_time_limit, ks_max_time_limit, ks_decay_factor, ks_feasibility_emphasis, true);
		}

		if (solve_alns)
		{
			auto initial_solution = InitalSolutionGenerator::GenerateInitialSolution(inst, true);
			ALNS alns(inst, initial_solution, alns_pool_size);

			delete initial_solution;
			initial_solution = nullptr;

			std::cout << MetaHeuristicSolution::GenerateALNSFileName(alns_num_iterations, alns_pool_size, metaheuristic_multi_threading) << "_seed_" << seed << std::endl;

			alns.Run(alns_num_iterations, metaheuristic_multi_threading);
			std::cout << "profits: " << (alns.best_solution())->profits_sum_ << std::endl;
			initial_sol = new MetaHeuristicSolution(alns.best_solution());
		}

		if (solve_simulated_annealing)
		{
			auto initial_solution = InitalSolutionGenerator::GenerateInitialSolution(inst, true);
			SimulatedAnnealing sa(inst, initial_solution);

			delete initial_solution;
			initial_solution = nullptr;

			std::cout << MetaHeuristicSolution::GenerateSimulatedAnnealingFileName(sa_temperature_decay_rate, metaheuristic_multi_threading) << "_seed_" << seed << std::endl;

			sa.Run(sa_temperature_decay_rate, sa_sampling_size, sa_target_acceptance_probability, metaheuristic_multi_threading);
			std::cout << "profits: " << (sa.best_solution())->profits_sum_ << std::endl;
			initial_sol = new MetaHeuristicSolution(sa.best_solution());
		}

		// in case no warm start is given, default one solving the integer problem considering only the mandatory vertices.
		if (initial_sol == nullptr)
		{
			std::cout << "generating DEFAULT initial sol" << std::endl;
			initial_sol = InitalSolutionGenerator::GenerateInitialSolution(inst, true);
		}
	}

	if (solve_compact)
	{
		auto algo = Solution<double>::GenerateAlgorithmName(solve_relaxed, solve_cb, solve_benders, (*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT], (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT], baseline, capacity_based,
															combine_feas_opt_cuts, separate_benders_cuts_relaxation, solve_generic_callback);
		std::cout << algo << std::endl;
		if (baseline)
		{
			std::list<UserCutGeneral *> *root_cuts = nullptr;

			// compute relaxation just to get limits at root node.
			bool use_valid_inequalities = false;
			bool find_root_cuts = false;
			CompactBaseline(inst, R0, Rn, time_limit, true, use_valid_inequalities, find_root_cuts, nullptr, initial_sol, force_use_all_vehicles, export_model, root_cuts, sol);

			if (!solve_relaxed)
				CompactBaseline(inst, R0, Rn, time_limit, false, use_valid_inequalities, find_root_cuts, nullptr, initial_sol, force_use_all_vehicles, export_model, root_cuts, sol);

			sol.write_to_file(algo, dir_solutions, instance_name);
		}

		if (capacity_based)
		{
			std::list<UserCutGeneral *> *root_cuts = nullptr;

			// compute relaxation just to get limits at root node.
			bool use_valid_inequalities = false;
			bool find_root_cuts = false;
			CompactSingleCommodity(inst, R0, Rn, time_limit, true, use_valid_inequalities, find_root_cuts, nullptr, initial_sol, force_use_all_vehicles, export_model, root_cuts, sol);

			if (!solve_relaxed)
				CompactSingleCommodity(inst, R0, Rn, time_limit, false, use_valid_inequalities, find_root_cuts, nullptr, initial_sol, force_use_all_vehicles, export_model, root_cuts, sol);

			sol.write_to_file(algo, dir_solutions, instance_name);
		}
	}

	if (solve_cb && !solve_benders)
	{
		auto root_cuts = new std::list<UserCutGeneral *>();
		auto algo = Solution<double>::GenerateAlgorithmName(solve_relaxed, solve_cb, solve_benders, (*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT], (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT], baseline, capacity_based,
															combine_feas_opt_cuts, separate_benders_cuts_relaxation, solve_generic_callback);

		std::cout << algo << std::endl;
		if (baseline)
		{
			bool use_valid_inequalities = false;
			bool find_root_cuts = true;
			CompactBaseline(inst, R0, Rn, time_limit, true, use_valid_inequalities, find_root_cuts, nullptr, initial_sol, force_use_all_vehicles, export_model, root_cuts, sol);

			if (time_limit != -1)
				time_limit = std::max(0.0, time_limit - sol.root_time_);
			sol.milp_time_ = sol.root_time_;

			use_valid_inequalities = true;
			find_root_cuts = false;
			if (!solve_relaxed)
				CompactBaseline(inst, R0, Rn, time_limit, false, use_valid_inequalities, find_root_cuts, root_cuts, initial_sol, force_use_all_vehicles, export_model, nullptr, sol);

			sol.write_to_file(algo, dir_solutions, instance_name);
		}

		if (capacity_based)
		{
			bool use_valid_inequalities = false;
			bool find_root_cuts = true;

			CompactSingleCommodity(inst, R0, Rn, time_limit, true, use_valid_inequalities, find_root_cuts, nullptr, initial_sol, force_use_all_vehicles, export_model, root_cuts, sol);

			// std::cout << root_cuts->size() << std::endl;
			// getchar();
			// getchar();

			if (time_limit != -1)
				time_limit = std::max(0.0, time_limit - sol.root_time_);
			sol.milp_time_ = sol.root_time_;

			use_valid_inequalities = true;
			find_root_cuts = false;
			if (!solve_relaxed)
				CompactSingleCommodity(inst, R0, Rn, time_limit, false, use_valid_inequalities, find_root_cuts, root_cuts, initial_sol, force_use_all_vehicles, export_model, nullptr, sol);

			sol.write_to_file(algo, dir_solutions, instance_name);
		}

		DeleteCuts(root_cuts);
		if (root_cuts != nullptr)
		{
			delete root_cuts;
			root_cuts = nullptr;
		}
	}

	if (solve_benders)
	{
		auto algo = Solution<double>::GenerateAlgorithmName(solve_relaxed, solve_cb, solve_benders, (*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT], (*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT], baseline, capacity_based,
															combine_feas_opt_cuts, separate_benders_cuts_relaxation, solve_generic_callback);
		std::cout << algo << std::endl;
		auto root_cuts = new std::list<UserCutGeneral *>();

		Formulation formulation;

		if (baseline)
		{
			formulation = Formulation::baseline;
		}
		if (capacity_based)
		{
			formulation = Formulation::single_commodity;
		}

		bool use_valid_inequalities = false;
		bool find_root_cuts = false;
		if (solve_cb)
		{
			use_valid_inequalities = false;
			find_root_cuts = true;

			Benders(inst, formulation, R0, Rn, time_limit, solve_generic_callback, combine_feas_opt_cuts, separate_benders_cuts_relaxation, true, use_valid_inequalities, find_root_cuts, nullptr, initial_sol, force_use_all_vehicles, export_model, root_cuts, sol);

			if (time_limit != -1)
				time_limit = std::max(0.0, time_limit - sol.root_time_);
			sol.milp_time_ = sol.root_time_;

			use_valid_inequalities = true;
			find_root_cuts = false;
		}
		else if (solve_relaxed)
		{
			Benders(inst, formulation, R0, Rn, time_limit, solve_generic_callback, combine_feas_opt_cuts, separate_benders_cuts_relaxation, true, use_valid_inequalities, find_root_cuts, nullptr, initial_sol, force_use_all_vehicles, export_model, root_cuts, sol);
		}

		if (!solve_relaxed)
			Benders(inst, formulation, R0, Rn, time_limit, solve_generic_callback, combine_feas_opt_cuts, separate_benders_cuts_relaxation, false, use_valid_inequalities, find_root_cuts, root_cuts, initial_sol, force_use_all_vehicles, export_model, nullptr, sol);

		// std::cout << "# opt cuts: " << sol.num_benders_opt_cuts_ << std::endl;
		// std::cout << "# feas cuts: " << sol.num_benders_feas_cuts_ << std::endl;
		// std::cout << "num cuts: " << sol.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << sol.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

		sol.write_to_file(algo, dir_solutions, instance_name);

		DeleteCuts(root_cuts);
		if (root_cuts != nullptr)
		{
			delete root_cuts;
			root_cuts = nullptr;
		}
	}

	if (!compute_initial_solution_heuristic && solve_kernel_search)
	{
		Formulation formulation;
		if (baseline)
		{
			formulation = Formulation::baseline;
		}
		if (capacity_based)
		{
			formulation = Formulation::single_commodity;
		}

		KernelSearch ks(inst);
		std::cout << KSHeuristicSolution::GenerateFileName(formulation, ks_max_size_bucket, ks_min_time_limit, ks_max_time_limit, ks_decay_factor, ks_feasibility_emphasis) << std::endl;
		const auto *kernel_search_sol = ks.Run(formulation, ks_max_size_bucket, ks_min_time_limit, ks_max_time_limit, ks_decay_factor, ks_feasibility_emphasis, true);
		kernel_search_sol->WriteToFile(inst, KSHeuristicSolution::GenerateFileName(formulation, ks_max_size_bucket, ks_min_time_limit, ks_max_time_limit, ks_decay_factor, ks_feasibility_emphasis), dir_solutions, instance_name);
		if (kernel_search_sol->is_feasible_)
			std::cout << " lb: " << kernel_search_sol->profits_sum_ << std::endl;
		if (kernel_search_sol->is_infeasible_)
			std::cout << " infeasible " << std::endl;
		std::cout << "time: " << kernel_search_sol->total_time_spent_ << std::endl;

		// if (kernel_search_sol->is_feasible_)
		// 	std::cout << -kernel_search_sol->profits_sum_ << std::endl;
		// else
		// 	std::cout << 0.0 << std::endl;

		delete kernel_search_sol;
		kernel_search_sol = nullptr;
	}

	if (!compute_initial_solution_heuristic && solve_alns)
	{
		auto initial_solution = InitalSolutionGenerator::GenerateInitialSolution(inst, true);
		ALNS alns(inst, initial_solution, alns_pool_size);

		delete initial_solution;
		initial_solution = nullptr;

		std::cout << MetaHeuristicSolution::GenerateALNSFileName(alns_num_iterations, alns_pool_size, metaheuristic_multi_threading) << "_seed_" << seed << std::endl;

		alns.Run(alns_num_iterations, metaheuristic_multi_threading);
		std::cout << "profits: " << (alns.best_solution())->profits_sum_ << std::endl;
		std::cout << "last improve iter: " << (alns.best_solution())->last_improve_iteration_ << std::endl;
		std::cout << "total iter: " << (alns.best_solution())->num_iterations_ << std::endl;
		std::cout << "time: " << (alns.best_solution())->total_time_spent_ << std::endl;
		(alns.best_solution())->WriteToFile(inst, MetaHeuristicSolution::GenerateALNSFileName(alns_num_iterations, alns_pool_size, metaheuristic_multi_threading) + "_seed_" + std::to_string(seed), dir_solutions, instance_name);
	}

	if (!compute_initial_solution_heuristic && solve_simulated_annealing)
	{
		auto initial_solution = InitalSolutionGenerator::GenerateInitialSolution(inst, true);
		SimulatedAnnealing sa(inst, initial_solution);

		delete initial_solution;
		initial_solution = nullptr;

		std::cout << MetaHeuristicSolution::GenerateSimulatedAnnealingFileName(sa_temperature_decay_rate, metaheuristic_multi_threading) << "_seed_" << seed << std::endl;

		sa.Run(sa_temperature_decay_rate, sa_sampling_size, sa_target_acceptance_probability, metaheuristic_multi_threading);
		// std::cout << *sa.best_solution() << std::endl;
		if (sa.best_solution()->is_infeasible_)
			std::cout << "INFEASIBLE" << std::endl;
		std::cout << "profits: " << (sa.best_solution())->profits_sum_ << std ::endl;
		std::cout << "last improve iter: " << (sa.best_solution())->last_improve_iteration_ << std::endl;
		std::cout << "total iter: " << (sa.best_solution())->num_iterations_ << std::endl;
		std::cout << "time: " << (sa.best_solution())->total_time_spent_ << std::endl;
		(sa.best_solution())->WriteToFile(inst, MetaHeuristicSolution::GenerateSimulatedAnnealingFileName(sa_temperature_decay_rate, metaheuristic_multi_threading) + "_seed_" + std::to_string(seed), dir_solutions, instance_name);
		// std::cout << *(sa.best_solution()) << std::endl;
	}

	if (solve_compact || solve_cb || solve_benders)
	{
		if (solve_relaxed)
		{
			std::cout << " LP: " << sol.lp_ << std::endl;
		}
		else
		{
			if (!sol.is_feasible_)
			{
				std::cout << " * infeasible" << std::endl;
			}
			else
			{
				sol.is_optimal_ ? std::cout << " optimal: " << sol.lb_ << std::endl
								: std::cout << " non optimal: [" << sol.lb_ << ", " << sol.ub_ << "]" << std::endl;
			}
		}
		std::cout << "num cuts: " << sol.num_cuts_added_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << "/" << sol.num_cuts_found_lp_[K_TYPE_CLIQUE_CONFLICT_CUT] << std::endl;

		if (solve_benders)
		{
			std::cout << "# opt cuts: " << sol.num_benders_opt_cuts_ << std::endl;
			std::cout << "# feas cuts: " << sol.num_benders_feas_cuts_ << std::endl;
		}
		std::cout << "root time: " << sol.root_time_ << std::endl;
		std::cout << "milp time: " << sol.milp_time_ << std::endl;
	}

	delete[] R0;
	R0 = nullptr;

	delete[] Rn;
	Rn = nullptr;

	if (initial_sol)
	{
		delete initial_sol;
		initial_sol = nullptr;
	}
}

int main(int argc, char *argv[])
{
	try
	{
		if (K_GETOPT)
		{
			ParseArgumentsAndRun(argc, argv);
			DeleteTimer();
			DeleteCallbackSelection();
		}
	}
	catch (const std::runtime_error &re)
	{
		std::cout << "Runtime error: " << re.what() << std::endl;
	}
	catch (const std::exception &ex)
	{
		std::cout << "Error occurred: " << ex.what() << std::endl;
	}
	catch (const int &error)
	{
		std::cout << "Error occurred: " << error << std::endl;
	}
	catch (IloException &e)
	{
		std::cout << "Concert Exception: " << e << std::endl;
	}
	catch (const char *e)
	{
		std::cout << e << std::endl;
	}
	catch (const std::string &e)
	{
		std::cout << e << std::endl;
	}
	catch (...)
	{
		std::cout << "Unknown failure occurred. Possible memory corruption" << std::endl;
	}

	return 0;
}

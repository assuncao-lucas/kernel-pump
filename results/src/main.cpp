#include <iostream>
#include <string>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <fmt/format.h>
#include <utils/floats.h>
#include "kernelpump/mipmodel.h"
#include "results/src/csv_reader.hpp"

#ifdef HAS_CPLEX
#include "kernelpump/cpxmodel.h"
#endif
#ifdef HAS_XPRESS
#include "kernelpump/xprsmodel.h"
#endif
#ifdef HAS_SCIP
#include "kernelpump/scipmodel.h"
#endif
#if defined(HAS_SCIP) && defined(HAS_ORTOOLS)
#include "kernelpump/pdlpmodel.h"
#endif

std::unordered_map<std::string, std::string> map_config_latex_name{
    {"config1", "FP$^*$"},
    {"config3", "FP$^{-}$"},
    {"config20", "FP$^{+}$"},
    {"config16", "KP$^*$"},
    {"config17", "KP$^{-}$"},
    {"config21", "KP$^{+}$"},
    {"config23", "KP$^{+/-}$"},
    {"config24", "KP$^{0}$"},
    {"config19", "CPLEX$_{std}$"},
    {"config18", "CPLEX$_{feas}$"},
};

std::unordered_set<std::string> map_algo_is_kernel_pump{
    "config16", "config17", "config21", "config23"};

std::unordered_set<std::string> map_algo_is_cplex{
    "config19", "config18"};

std::tuple<bool, double> is_number(const std::string &s)
{
    char *end = nullptr;
    double val = strtod(s.c_str(), &end);
    return {end != s.c_str() && *end == '\0' && val != HUGE_VAL, val};
}

void AddFiles(std::string directory, std::vector<std::pair<std::string, std::string>> &files)
{
    std::ifstream file(directory.c_str());
    std::string curr_file, curr_format_termination;

    if (file.is_open())
    {
        while (std::getline(file, curr_file))
        {
            remove_if(curr_file.begin(), curr_file.end(), isspace);
            if (!curr_file.empty())
            {
                size_t end_pos = curr_file.find_last_of(".");

                if (end_pos != std::string::npos)
                {
                    curr_format_termination = curr_file.substr(end_pos, std::string::npos);
                    curr_file = curr_file.substr(0, end_pos);
                }

                files.push_back({curr_file, curr_format_termination});
                // std::cout << curr_file << " " << curr_format_termination << std::endl;
            }
        }
    }
    else
    {
        std::cout << "Invalid directory" << std::endl;
        throw 1;
    }
}

void AddFilesFromDirectory(std::string directory, std::vector<std::pair<std::string, std::string>> &files, bool add_dir)
{
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(directory.c_str())) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL)
        {
            std::string curr_file(ent->d_name);
            std::string curr_format_termination;
            if ((curr_file != ".") && (curr_file != "..") && ((curr_file.size() > 0) && (curr_file[curr_file.size() - 1] != '~')))
            {
                size_t end_pos = curr_file.find(".");
                if (end_pos != std::string::npos)
                {
                    curr_format_termination = curr_file.substr(end_pos, std::string::npos);
                    curr_file = curr_file.substr(0, end_pos);
                }

                std::string inst_path;
                if (add_dir)
                    inst_path = directory;
                inst_path.append(curr_file);
                files.push_back({inst_path, curr_format_termination});
                // std::cout << curr_file << " " << curr_format_termination << std::endl;
            }
        }
        closedir(dir);
    }
    else
    {
        std::cout << "Invalid directory" << std::endl;
        throw 1;
    }
    // std::sort(instances.begin(),instances.end());
}

enum class PerformanceMeasureType
{
    success,
    time
};

struct PerformanceTuple
{
    int wins;
    int ties;
    int losses;
};

void GeneratePerformanceMatrix(std::string instances_folder, std::string inst_list_file, std::string solutions_folder, double time_limit, PerformanceMeasureType type)
{
    std::string curr_file;
    std::vector<std::pair<std::string, std::string>> instances;
    std::vector<std::string> configs{
        "config1", "config3", "config20", "config16", "config17", "config21", "config23", "config19", "config18"};
    std::vector<std::string>
        seeds{"1", "2", "3", "4", "5"};
    int num_seeds = seeds.size();

    // AddFilesFromDirectory(instances_folder, instances, false);
    AddFiles(inst_list_file, instances);

    std::sort(instances.begin(), instances.end(), [](const auto &a, const auto &b)
              { return a.first < b.first; });
    // std::sort(configs.begin(), configs.end());

    // for (const auto &[key, value] : instances_bounds)
    //     std::cout << key << " " << value << std::endl;

    std::fstream output;
    std::string output_name;
    output_name = "..//tables//latex//performance_matrix.txt";
    output.open(output_name.c_str(), std::fstream::out);

    if (!output.is_open())
    {
        std::cout << "Could not open file " << output_name << std::endl;
        throw 1;
    }

    // std::cout << output_name << std::endl;

    output << std::setprecision(4) << std::fixed;

    int total_num_instances = instances.size();

    std::vector<std::vector<PerformanceTuple>> matrix(configs.size());

    output << " &";
    for (size_t config_num = 0; config_num < configs.size(); ++config_num)
    {
        matrix[config_num].resize(configs.size(), {0, 0, 0});

        output << map_config_latex_name[configs[config_num]];
        if (config_num < configs.size() - 1)
            output << " & ";
        else
            output << " \\\\" << std::endl;
    }

    for (const auto instance : instances)
    {
        std::vector<std::pair<double, double>> instance_results_per_config(configs.size(), {0, 0});
        for (size_t config_num = 0; config_num < configs.size(); ++config_num)
        {
            double avg_time = 0.0, avg_iter = 0.0, avg_success = 0.0;
            for (size_t seed_num = 0; seed_num < num_seeds; ++seed_num)
            {
                curr_file = solutions_folder + "//";
                curr_file.append("s_");
                curr_file.append(configs[config_num]);
                curr_file.append("_");
                curr_file.append(instance.first);
                curr_file.append("_");
                curr_file.append(seeds[seed_num]);
                curr_file.append(".sol");

                // std::cout << curr_file << std::endl;

                std::fstream input;
                input.open(curr_file.c_str(), std::fstream::in);

                if (!input.is_open())
                {
                    std::cout << "Could not open file " << curr_file << std::endl;
                    continue;
                }

                std::stringstream s_time, s_iter;
                std::string status;
                double iter = 0.0, total_time = 0.0;
                std::string line;

                getline(input, line);
                size_t pos = line.find_first_of(":");
                status = line.substr(pos + 2);

                std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
                status.erase(end_pos, status.end());

                getline(input, line);
                getline(input, line);
                pos = line.find_first_of(":");
                s_time << line.substr(pos + 2);
                s_time >> total_time;

                double tolerance = 10.0;
                if (greaterThan(total_time, time_limit + tolerance))
                {
                    std::cout << "* " << instance.first << " seed " << seed_num + 1 << " " << configs[config_num] << ", total time of " << total_time << " considered as " << time_limit << std::endl;
                    total_time = time_limit;
                }

                total_time = std::min(total_time, time_limit);

                getline(input, line);
                pos = line.find_first_of(":");
                s_iter << line.substr(pos + 2);
                s_iter >> iter;

                avg_time += total_time;
                avg_iter += iter;

                if (status == "FOUNDINTEGERFEASIBLE")
                {
                    avg_success += 1;
                }

                input.close();
            }

            avg_time /= num_seeds;
            avg_iter /= num_seeds;
            avg_success /= num_seeds;

            // if (type == PerformanceMeasureType::time)
            instance_results_per_config[config_num].second = avg_time;
            // else if (type == PerformanceMeasureType::success)
            instance_results_per_config[config_num].first = 1.0 - ceil(avg_success);
        }

        for (size_t i = 0; i < configs.size(); ++i)
        {
            for (size_t j = 0; j < configs.size(); ++j)
            {
                if (type == PerformanceMeasureType::success)
                {
                    if (lessThan(instance_results_per_config[i].first, instance_results_per_config[j].first))
                    {
                        matrix[i][j].wins += 1;
                        if (configs[i] == "config17" && configs[j] == "config19")
                            std::cout << "!!!!!!!!!!!!!!!!!!" << instance.first << std::endl;
                    }
                    else if (equal(instance_results_per_config[i].first, instance_results_per_config[j].first))
                        matrix[i][j].ties += 1;
                    else
                        matrix[i][j].losses += 1;
                }
                else if (type == PerformanceMeasureType::time)
                {
                    // // only consider the cases where both algorithms succeed in finding a feasible solution
                    // if (equal(instance_results_per_config[i].first, instance_results_per_config[j].first) && equal(instance_results_per_config[i].first, 0.0))
                    {
                        if (lessThan(instance_results_per_config[i].second, instance_results_per_config[j].second))
                            matrix[i][j].wins += 1;
                        else if (equal(instance_results_per_config[i].second, instance_results_per_config[j].second))
                            matrix[i][j].ties += 1;
                        else
                            matrix[i][j].losses += 1;
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < configs.size(); ++i)
    {
        output << map_config_latex_name[configs[i]] << " & ";

        for (size_t j = 0; j < configs.size(); ++j)
        {
            if (i == j)
                output << "{--}";
            else
                output << matrix[i][j].wins << "/" << matrix[i][j].ties << "/" << matrix[i][j].losses;

            if (j < configs.size() - 1)
                output << " & ";
            else
                output << " \\\\" << std::endl;
        }
    }

    output.close();
}

void GeneratePerformanceProfile(std::string instances_folder, std::string inst_list_file, std::string solutions_folder, double time_limit, PerformanceMeasureType type)
{
    std::string curr_file;
    std::vector<std::pair<std::string, std::string>> instances;
    std::vector<std::string> configs{
        "config1", "config16", "config17", "config21", "config23", "config18"};
    std::vector<std::string>
        seeds{"1", "2", "3", "4", "5"};
    int num_seeds = seeds.size();

    // AddFilesFromDirectory(instances_folder, instances, false);
    AddFiles(inst_list_file, instances);

    std::sort(instances.begin(), instances.end(), [](const auto &a, const auto &b)
              { return a.first < b.first; });
    // std::sort(configs.begin(), configs.end());

    // for (const auto &[key, value] : instances_bounds)
    //     std::cout << key << " " << value << std::endl;

    std::fstream output;
    std::string output_name;
    output_name = "..//tables//latex//performance_profile.txt";
    output.open(output_name.c_str(), std::fstream::out);

    if (!output.is_open())
    {
        std::cout << "Could not open file " << output_name << std::endl;
        throw 1;
    }

    // std::cout << output_name << std::endl;

    output << std::setprecision(4) << std::fixed;

    int total_num_instances = instances.size();

    for (size_t config_num = 0; config_num < configs.size(); ++config_num)
    {
        output << map_config_latex_name[configs[config_num]];
        if (config_num < configs.size() - 1)
            output << " ";
        else
            output << std::endl;
    }

    for (const auto instance : instances)
    {
        for (size_t config_num = 0; config_num < configs.size(); ++config_num)
        {
            double avg_time = 0.0, avg_iter = 0.0, avg_success = 0.0;
            for (size_t seed_num = 0; seed_num < num_seeds; ++seed_num)
            {
                curr_file = solutions_folder + "//";
                curr_file.append("s_");
                curr_file.append(configs[config_num]);
                curr_file.append("_");
                curr_file.append(instance.first);
                curr_file.append("_");
                curr_file.append(seeds[seed_num]);
                curr_file.append(".sol");

                // std::cout << curr_file << std::endl;

                std::fstream input;
                input.open(curr_file.c_str(), std::fstream::in);

                if (!input.is_open())
                {
                    std::cout << "Could not open file " << curr_file << std::endl;
                    continue;
                }

                std::stringstream s_time, s_iter;
                std::string status;
                double iter = 0.0, total_time = 0.0;
                std::string line;

                getline(input, line);
                size_t pos = line.find_first_of(":");
                status = line.substr(pos + 2);

                std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
                status.erase(end_pos, status.end());

                getline(input, line);
                getline(input, line);
                pos = line.find_first_of(":");
                s_time << line.substr(pos + 2);
                s_time >> total_time;

                double tolerance = 10.0;
                if (greaterThan(total_time, time_limit + tolerance))
                {
                    std::cout << "* " << instance.first << " seed " << seed_num + 1 << " " << configs[config_num] << ", total time of " << total_time << " considered as " << time_limit << std::endl;
                    total_time = time_limit;
                }

                total_time = std::min(total_time, time_limit);

                getline(input, line);
                pos = line.find_first_of(":");
                s_iter << line.substr(pos + 2);
                s_iter >> iter;

                avg_time += total_time;
                avg_iter += iter;

                if (status == "FOUNDINTEGERFEASIBLE")
                {
                    avg_success += 1;
                }

                input.close();
            }

            avg_time /= num_seeds;
            avg_iter /= num_seeds;
            avg_success /= num_seeds;

            if (type == PerformanceMeasureType::time)
                output << avg_time;

            if (config_num < configs.size() - 1)
                output << " ";
            else
                output << std::endl;
        }
    }

    output.close();
}

double StDev(const std::vector<double> &gaps)
{
    double stdev = 0.0;
    double b_avg_gap = 0.0;

    for (size_t i = 0; i < gaps.size(); i++)
    {
        b_avg_gap += gaps.at(i);
    }

    b_avg_gap /= (1.0 * gaps.size());

    for (size_t i = 0; i < gaps.size(); i++)
    {
        stdev += pow(gaps.at(i) - b_avg_gap, 2.0);
    }

    stdev = stdev / (gaps.size() * 1.0 - 1.0);
    return pow(stdev, 0.5);
}

void GenerateAlgorithmsCSVAndLatexTable(std::string instances_folder, std::string inst_list_file, std::string solutions_folder, std::string best_known_bounds_csv, double time_limit)
{
    std::string curr_file;
    std::string solver = "cpx";
    std::vector<std::pair<std::string, std::string>> instances;
    std::vector<std::string> configs{
        "config1", "config3", "config20", "config16", "config17", "config21", "config23", "config19", "config18"};
    std::vector<std::string>
        seeds{"1", "2", "3", "4", "5"};
    int num_seeds = seeds.size();

    MIPModelPtr model;
#ifdef HAS_CPLEX
    if (solver == "cpx")
        model = MIPModelPtr(new CPXModel());
#else
    if (solver == "cpx")
        throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif
#ifdef HAS_XPRESS
    if (solver == "xprs")
        model = MIPModelPtr(new XPRSModel());
#else
    if (solver == "xprs")
        throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif
#ifdef HAS_SCIP
    if (solver == "scip")
        model = MIPModelPtr(new SCIPModel());
#else
    if (solver == "scip")
        throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif
#if defined(HAS_SCIP) && defined(HAS_ORTOOLS)
    if (solver == "pdlp")
        model = MIPModelPtr(new PDLPModel());
#else
    if (solver == "pdlp")
        throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif

    if (!model)
        throw std::runtime_error("No solver available");

    // AddFilesFromDirectory(instances_folder, instances, false);
    AddFiles(inst_list_file, instances);

    std::sort(instances.begin(), instances.end(), [](const auto &a, const auto &b)
              { return a.first < b.first; });
    // std::sort(configs.begin(), configs.end());

    std::ifstream file(best_known_bounds_csv);
    std::unordered_map<std::string, double> instances_bounds;

    for (auto &row : CSVRange(file))
    {
        std::string instance_name(row[0]);
        std::string bound_str(row[10]);

        instance_name.erase(std::remove(instance_name.begin(), instance_name.end(), ' '), instance_name.end());
        bound_str.erase(std::remove(bound_str.begin(), bound_str.end(), ' '), bound_str.end());
        bound_str.erase(std::remove(bound_str.begin(), bound_str.end(), '*'), bound_str.end());

        auto [isValidNumber, value] = is_number(bound_str);
        if (isValidNumber)
            instances_bounds[instance_name] = value;
    }

    // for (const auto &[key, value] : instances_bounds)
    //     std::cout << key << " " << value << std::endl;

    std::fstream output, output2, output3;
    std::string output_name, output_name2, output_name3;
    output_name = "..//tables//latex//table_algorithms.csv";
    output_name2 = "..//tables//latex//table_performance.txt";
    output_name3 = "..//tables//latex//table_convergence.txt";
    output.open(output_name.c_str(), std::fstream::out);
    output2.open(output_name2.c_str(), std::fstream::out);
    output3.open(output_name3.c_str(), std::fstream::out);

    if (!output.is_open())
    {
        std::cout << "Could not open file " << output_name << std::endl;
        throw 1;
    }

    if (!output2.is_open())
    {
        std::cout << "Could not open file " << output_name2 << std::endl;
        throw 1;
    }

    if (!output3.is_open())
    {
        std::cout << "Could not open file " << output_name2 << std::endl;
        throw 1;
    }
    // std::cout << output_name << std::endl;

    output << std::setprecision(6) << std::fixed;
    output2 << std::setprecision(2) << std::fixed;
    output3 << std::setprecision(2) << std::fixed;

    std::vector<double> total_time_per_config(configs.size(), 0.0);
    std::vector<double> total_iter_per_config(configs.size(), 0.0);
    std::vector<double> total_proj_gap_per_config(configs.size(), 0.0);
    std::vector<double> total_actual_gap_per_config(configs.size(), 0.0);
    std::vector<double> total_obj_gap_per_config(configs.size(), 0.0);
    std::vector<double> total_avg_success_per_config(configs.size(), 0.0);
    std::vector<double> total_num_success_per_config(configs.size(), 0.0);
    std::vector<double> total_num_frac_per_config(configs.size(), 0.0);
    std::vector<double> total_percentage_visited_buckets_per_config(configs.size(), 0.0);
    std::vector<double> total_percentage_first_feasible_bucket_per_config(configs.size(), 0.0);
    std::vector<double> total_percentage_added_vars_per_config(configs.size(), 0.0);
    std::vector<double> total_percentage_active_vars_per_config(configs.size(), 0.0);
    std::vector<double> num_instances_discarded_from_obj_gap_computation_per_config(configs.size(), 0);

    std::vector<std::vector<double>> percentage_visited_buckets_per_config(configs.size());
    std::vector<std::vector<double>> percentage_added_vars_per_config(configs.size());
    std::vector<std::vector<double>> percentage_active_vars_per_config(configs.size());

    int total_num_instances = instances.size();
    int num_integer_and_binary_vars = 0;
    int num_binary_vars = 0;

    output << ", ";
    for (int i = 0; i < configs.size(); ++i)
        output << ", " << configs[i] << " , , , , , , , , , , , , ";
    // output << "\\\\" << std::endl;
    output << std::endl;

    output << "instance ";
    for (int i = 0; i < configs.size(); ++i)
        output << ", , total success, avg success , iter , visited buckets (%), first feasible bucket (%) , actual int gap (%) , proj int gap (%), num frac (%), obj gap (%), time (s) , added bin vars (%), active bin vars (%)";
    // output << "\\\\" << std::endl;
    output << std::endl;

    for (const auto instance : instances)
    {
        double inst_best_bound = instances_bounds[instance.first];
        // std::cout << instance.first + instance.second << std::endl;
        output << instance.first;
        model->readModel(instances_folder + "//" + instance.first + instance.second);
        // double objsen = (model->objSense() == ObjSense::MIN) ? 1.0 : -1.0;
        num_integer_and_binary_vars = model->numIntegerAndBinaryCols();
        num_binary_vars = model->numBinaryCols();
        // std::cout << "num int vars: " << num_integer_and_binary_vars << std::endl;
        for (size_t config_num = 0; config_num < configs.size(); ++config_num)
        {
            double avg_time = 0.0, avg_iter = 0.0, avg_proj_gap = 0.0, avg_actual_gap = 0.0, avg_success = 0.0, avg_num_frac = 0.0, avg_obj_gap = 0.0;
            double avg_percentage_visited_buckets = 0.0, avg_percentage_first_feasible_bucket = 0.0;
            double avg_percentage_added_vars = 0.0, avg_percentage_active_vars = 0.0;
            int num_exec_discarded_from_obj_gap_computation = 0;
            double obj_gap = 0.0;
            for (size_t seed_num = 0; seed_num < num_seeds; ++seed_num)
            {
                curr_file = solutions_folder + "//";
                curr_file.append("s_");
                curr_file.append(configs[config_num]);
                curr_file.append("_");
                curr_file.append(instance.first);
                curr_file.append("_");
                curr_file.append(seeds[seed_num]);
                curr_file.append(".sol");

                // std::cout << curr_file << std::endl;

                std::fstream input;
                input.open(curr_file.c_str(), std::fstream::in);

                if (!input.is_open())
                {
                    std::cout << "Could not open file " << curr_file << std::endl;
                    continue;
                }

                std::stringstream s_time, s_iter, s_num_buckets, s_last_visited_bucket, s_first_feasible_bucket, s_proj_gap, s_actual_gap, s_num_frac, s_obj_value;
                std::stringstream s_num_added_bin_vars, s_num_active_bin_vars;
                std::string status;
                double iter = 0.0, total_time = 0.0, time_build = 0.0, proj_gap = 0.0, actual_gap = 0.0, num_frac = 0.0, obj_value;
                double num_buckets = 0.0, last_visited_bucket = 0.0, first_feasible_bucket = 0.0;
                double num_added_bin_vars = 0.0, num_active_bin_vars = 0.0;
                std::string line;

                getline(input, line);
                size_t pos = line.find_first_of(":");
                status = line.substr(pos + 2);

                std::string::iterator end_pos = std::remove(status.begin(), status.end(), ' ');
                status.erase(end_pos, status.end());

                getline(input, line);
                pos = line.find_first_of(":");
                s_time << line.substr(pos + 2);
                s_time >> time_build;
                s_time.clear();

                getline(input, line);
                pos = line.find_first_of(":");
                s_time << line.substr(pos + 2);
                s_time >> total_time;

                double tolerance = 10.0;
                if (greaterThan(total_time, time_limit + tolerance))
                {
                    std::cout << "* " << instance.first << " seed " << seed_num + 1 << " " << configs[config_num] << ", total time of " << total_time << " considered as " << time_limit << std::endl;
                    total_time = time_limit;
                }

                total_time = std::min(total_time, time_limit);

                getline(input, line);
                pos = line.find_first_of(":");
                s_iter << line.substr(pos + 2);
                s_iter >> iter;

                getline(input, line);
                pos = line.find_first_of(":");
                s_num_buckets << line.substr(pos + 2);
                s_num_buckets >> num_buckets;

                getline(input, line);
                pos = line.find_first_of(":");
                s_last_visited_bucket << line.substr(pos + 2);
                s_last_visited_bucket >> last_visited_bucket;

                getline(input, line);
                pos = line.find_first_of(":");
                s_first_feasible_bucket << line.substr(pos + 2);
                s_first_feasible_bucket >> first_feasible_bucket;

                for (int i = 0; i < 2; ++i)
                    getline(input, line);

                pos = line.find_first_of(":");
                s_obj_value << line.substr(pos + 2);
                s_obj_value >> obj_value;

                getline(input, line);
                pos = line.find_first_of(":");
                s_actual_gap << line.substr(pos + 2);
                s_actual_gap >> actual_gap;

                getline(input, line);
                pos = line.find_first_of(":");
                s_proj_gap << line.substr(pos + 2);
                s_proj_gap >> proj_gap;

                if (equal(proj_gap, INFBOUND)) // handle case when projection bound is at its default value (when no pump iter was completed).
                    proj_gap = num_integer_and_binary_vars;

                getline(input, line);
                pos = line.find_first_of(":");
                s_num_frac << line.substr(pos + 2);
                s_num_frac >> num_frac;

                if (getline(input, line))
                {
                    if (!line.empty())
                    {
                        pos = line.find_first_of(":");
                        s_num_added_bin_vars << line.substr(pos + 2);
                        s_num_added_bin_vars >> num_added_bin_vars;
                    }
                }

                if (getline(input, line))
                {
                    if (!line.empty())
                    {
                        pos = line.find_first_of(":");
                        s_num_active_bin_vars << line.substr(pos + 2);
                        s_num_active_bin_vars >> num_active_bin_vars;
                    }
                }

                avg_time += total_time;
                avg_iter += iter;
                avg_proj_gap += 100.0 * (proj_gap / num_integer_and_binary_vars);
                avg_actual_gap += 100.0 * (actual_gap / num_integer_and_binary_vars);
                avg_num_frac += 100.0 * (num_frac / num_integer_and_binary_vars);
                avg_percentage_first_feasible_bucket += 100.0 * ((last_visited_bucket == -1 || first_feasible_bucket == -1 || first_feasible_bucket == num_buckets) ? 1.0 : (first_feasible_bucket + 1) / (num_buckets + 1));

                obj_gap = 0.0;
                // std::cout << "obj value " << obj_value << " best bound: " << inst_best_bound << std::endl;
                if (status == "FOUNDINTEGERFEASIBLE")
                {
                    avg_success += 1;
                    // compute obj value gap
                    if (equal(obj_value, inst_best_bound))
                    {
                        obj_gap = 0.0;
                    }
                    else if (lessThan(inst_best_bound * obj_value, 0.0))
                    {
                        obj_gap = 100.0;
                    }
                    else
                    {
                        obj_gap = (fabs(inst_best_bound - obj_value) / std::max(fabs(inst_best_bound), fabs(obj_value))) * 100.0;
                    }

                    avg_percentage_visited_buckets += 100.0 * ((last_visited_bucket == -1 || last_visited_bucket == num_buckets) ? 1.0 : (last_visited_bucket + 1) / (num_buckets + 1));
                    avg_percentage_active_vars += (num_binary_vars == 0) ? 100.0 : 100.0 * num_active_bin_vars / num_binary_vars;
                    avg_percentage_added_vars += (num_binary_vars == 0) ? 100.0 : 100.0 * num_added_bin_vars / num_binary_vars;
                }
                else
                {
                    ++num_exec_discarded_from_obj_gap_computation;
                    // std::cout << " *** discarded from obj gap computation due to failure in finding integer solution" << std::endl;
                }
                avg_obj_gap += obj_gap;
                // std::cout << "config " << config_num << " gap: " << obj_gap << std::endl;

                input.close();
            }

            if (num_exec_discarded_from_obj_gap_computation < num_seeds)
            {
                avg_obj_gap /= (num_seeds - num_exec_discarded_from_obj_gap_computation);
                avg_percentage_active_vars /= (num_seeds - num_exec_discarded_from_obj_gap_computation);
                avg_percentage_added_vars /= (num_seeds - num_exec_discarded_from_obj_gap_computation);
                avg_percentage_visited_buckets /= (num_seeds - num_exec_discarded_from_obj_gap_computation);

                total_obj_gap_per_config[config_num] += avg_obj_gap;
                total_percentage_added_vars_per_config[config_num] += avg_percentage_added_vars;
                total_percentage_active_vars_per_config[config_num] += avg_percentage_active_vars;
                total_percentage_visited_buckets_per_config[config_num] += avg_percentage_visited_buckets;

                percentage_added_vars_per_config[config_num].push_back(avg_percentage_added_vars);
                percentage_active_vars_per_config[config_num].push_back(avg_percentage_active_vars);
                percentage_visited_buckets_per_config[config_num].push_back(avg_percentage_visited_buckets);
            }
            else // if execution of that instance discarded for all seeds, do not take it into account in the total computation.
            {
                (num_instances_discarded_from_obj_gap_computation_per_config[config_num]) += 1;
                avg_obj_gap = -1.0;
                avg_percentage_active_vars = -1.0;
                avg_percentage_added_vars = -1.0;
                avg_percentage_visited_buckets = -1.0;
            }
            avg_time /= num_seeds;
            avg_iter /= num_seeds;
            avg_proj_gap /= num_seeds;
            avg_actual_gap /= num_seeds;
            avg_success /= num_seeds;
            avg_num_frac /= num_seeds;
            avg_percentage_first_feasible_bucket /= num_seeds;

            output << " , , " << std::ceil(avg_success) << " , " << avg_success << " , " << avg_iter << " , " << avg_percentage_visited_buckets << " , " << avg_percentage_first_feasible_bucket << " , " << avg_actual_gap << " , " << avg_proj_gap << " , " << avg_num_frac << " , " << avg_obj_gap << " , " << avg_time << " , " << avg_percentage_added_vars << " , " << avg_percentage_active_vars;

            total_avg_success_per_config[config_num] += avg_success;
            total_num_success_per_config[config_num] += std::ceil(avg_success);
            total_iter_per_config[config_num] += avg_iter;
            total_proj_gap_per_config[config_num] += avg_proj_gap;
            total_actual_gap_per_config[config_num] += avg_actual_gap;
            total_time_per_config[config_num] += avg_time;
            total_num_frac_per_config[config_num] += avg_num_frac;
            total_percentage_first_feasible_bucket_per_config[config_num] += avg_percentage_first_feasible_bucket;
        }
        // output << "\\\\" << std::endl;
        output << std::endl;
    }

    output << "Total";
    for (size_t config_num = 0; config_num < configs.size(); ++config_num)
    {
        if (num_instances_discarded_from_obj_gap_computation_per_config[config_num] >= total_num_instances)
            output << " , , " << (int)(total_num_success_per_config[config_num]) << " , " << total_avg_success_per_config[config_num] / total_num_instances << " , " << total_iter_per_config[config_num] / total_num_instances << " , - , " << total_percentage_first_feasible_bucket_per_config[config_num] / total_num_instances << " , " << total_actual_gap_per_config[config_num] / total_num_instances << " , " << total_proj_gap_per_config[config_num] / total_num_instances << " , " << total_num_frac_per_config[config_num] / total_num_instances << " , - , " << total_time_per_config[config_num] / total_num_instances << " , - , - ";
        else
            output << " , , " << (int)(total_num_success_per_config[config_num])
                   << " , " << total_avg_success_per_config[config_num] / total_num_instances << " , " << total_iter_per_config[config_num] / total_num_instances << " , "
                   << total_percentage_visited_buckets_per_config[config_num] / (total_num_instances - num_instances_discarded_from_obj_gap_computation_per_config[config_num]) << " , " << total_percentage_first_feasible_bucket_per_config[config_num] / total_num_instances
                   << " , " << total_actual_gap_per_config[config_num] / total_num_instances << " , " << total_proj_gap_per_config[config_num] / total_num_instances << " , "
                   << total_num_frac_per_config[config_num] / total_num_instances << " , " << total_obj_gap_per_config[config_num] / (total_num_instances - num_instances_discarded_from_obj_gap_computation_per_config[config_num])
                   << " , " << total_time_per_config[config_num] / total_num_instances
                   << " , " << total_percentage_added_vars_per_config[config_num] / (total_num_instances - num_instances_discarded_from_obj_gap_computation_per_config[config_num])
                   << " , " << total_percentage_active_vars_per_config[config_num] / (total_num_instances - num_instances_discarded_from_obj_gap_computation_per_config[config_num]);

        output2 << map_config_latex_name[configs[config_num]] << " && "
                << (int)(total_num_success_per_config[config_num]) << "/" << total_num_instances
                << " & " << total_time_per_config[config_num] / total_num_instances
                << " & ";
        (map_algo_is_cplex.find(configs[config_num]) == map_algo_is_cplex.end()) ? output2 << total_num_frac_per_config[config_num] / total_num_instances : output2 << " {--} ";

        output2 << " \\\\" << std::endl;

        if (map_algo_is_kernel_pump.find(configs[config_num]) != map_algo_is_kernel_pump.end())
        {
            output3 << map_config_latex_name[configs[config_num]] << " && "
                    << total_percentage_visited_buckets_per_config[config_num] / (total_num_instances - num_instances_discarded_from_obj_gap_computation_per_config[config_num])
                    << " & " << StDev(percentage_visited_buckets_per_config[config_num])
                    << " & " << total_percentage_added_vars_per_config[config_num] / (total_num_instances - num_instances_discarded_from_obj_gap_computation_per_config[config_num])
                    << " & " << StDev(percentage_added_vars_per_config[config_num])
                    << " & " << total_percentage_active_vars_per_config[config_num] / (total_num_instances - num_instances_discarded_from_obj_gap_computation_per_config[config_num])
                    << " & " << StDev(percentage_active_vars_per_config[config_num])
                    << " \\\\" << std::endl;
        }
    }

    output << std::endl;

    output.close();
    output2.close();
    output3.close();
}

void ComputeDifferenceLists(std::string list1_file, std::string list2_file)
{
    std::vector<std::pair<std::string, std::string>> full_list1, full_list2;
    std::vector<std::string> list1, list2;
    AddFiles(list1_file, full_list1);
    AddFiles(list2_file, full_list2);

    for (const auto &element : full_list1)
        list1.push_back(element.first + element.second);

    for (const auto &element : full_list2)
        list2.push_back(element.first + element.second);

    std::sort(list1.begin(), list1.end());
    std::sort(list2.begin(), list2.end());

    std::vector<std::string> difference;

    // Step 5: Use set_difference
    set_difference(list1.begin(), list1.end(), list2.begin(),
                   list2.end(), back_inserter(difference));

    for (const auto &element : difference)
        std::cout << element << std::endl;
}

int main()
{
    // ComputeDifferenceLists("/home/lucas/Documents/Research/kernel-pump/instances/integer_problems_collection.txt", "/home/lucas/Documents/Research/kernel-pump/instances/integer_problems_benchmark.txt");
    std::string inst_folder = "/home/lucas/Downloads/instances/benchmark";
    std::string inst_list_file = "/home/lucas/Documents/Research/kernel-pump/instances/all_problems_benchmark.txt";
    std::string solutions_folder = "/home/lucas/Documents/Research/kernel-pump/solutions/miplib-new";
    std::string best_known_bounds_csv = "/home/lucas/Documents/Research/kernel-pump/results/bestKnownBounds.csv";

    // GenerateAlgorithmsCSVAndLatexTable(inst_folder, inst_list_file, solutions_folder, best_known_bounds_csv, 3600.0);

    // GeneratePerformanceProfile(inst_folder, inst_list_file, solutions_folder, 3600.0, PerformanceMeasureType::time);
    GeneratePerformanceMatrix(inst_folder, inst_list_file, solutions_folder, 3600.0, PerformanceMeasureType::success);
    return 0;
}
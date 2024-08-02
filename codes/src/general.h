#pragma once

#include <cmath>
#include <vector>
#include <iostream>
#include <utility>

const double K_PRECISION_COMPARE_DOUBLE = 0.0000001;

const std::string K_FILE_DELIMITER = "***********************************";

// feasibility pump parameters
const double K_PUMP_IMPROVEMENT_TOLERANCE = 0.1;
const int max_iter_stage1 = 1000;
const int max_iter_stage2 = 200;
const int max_stalls_stage1 = 10;
const int max_stalls_stage2 = 10;
const double alpha_decrease_rate = 0.9;
const double initial_alpha_stage1 = 0.0;
const double initial_alpha_stage2 = 0.0; // if zero, then feasibility pump, if 1.0, then objective feasibility pump.
const bool K_SOLVE_STAGE_1 = false;
const double perturbation_flip_percentage = 0.1;
const int K_FLIP_BASIS = 10;
const bool FIXED_FLIP_BASIS = true;
const double K_APLHA_DECREMENT_PRECISION = 0.005;
const bool K_FEASIBILITY_PUMP_ADD_CUTS = false;

// DEFAULT Kernel Search parameters
const bool K_KERNEL_SEARCH_ADD_CUTS = false;
const int K_KS_MAX_SIZE_BUCKET = 5;
const int K_KS_MAX_TIME_LIMIT = 900;
const int K_KS_MIN_TIME_LIMIT = 5;
const double K_KS_DECAY_FACTOR_TIME_LIMIT = 0.9;

double StDev(std::vector<double> &gaps, double b_avg_gap);
void split_file_path(std::string directory, std::string &folder, std::string &file_name);
bool double_equals(double a, double b, double epsilon = K_PRECISION_COMPARE_DOUBLE);
bool double_greater(double a, double b, double epsilon = K_PRECISION_COMPARE_DOUBLE);
bool double_less(double a, double b, double epsilon = K_PRECISION_COMPARE_DOUBLE);
double round_decimals(double value, int num_decimals);
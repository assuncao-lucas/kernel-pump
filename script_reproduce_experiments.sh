#!/bin/bash

# Default values
# Mandatory arguments (no defaults)
instances_dir=""
instances_list=""
configs_list=""
time_limit=""
num_seeds=""

# optional argument.
percentage_random_sample=100

base_solutions_dir="./solutions" 
solutions_sub_folder=$(date +%Y-%m-%d_%H:%M:%S)
solutions_dir="$base_solutions_dir/$solutions_sub_folder"
configs_dir="./settings/all"
usage_str="Usage: $0 --instances-dir PATH --instances-list FILE --configs-list FILE --time-limit FLOAT --num-seeds INT>0 [--percentage-random-sample INTEGER 0-100]"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --instances-dir)
            instances_dir="$2"
            shift 2
            ;;
        --instances-list)
            instances_list="$2"
            shift 2
            ;;
        --configs-list)
            configs_list="$2"
            shift 2
            ;;
        --percentage-random-sample)
            percentage_random_sample="$2"
            shift 2
            ;;
        --time-limit)
            time_limit="$2"
            shift 2
            ;;
        --num-seeds)
            num_seeds="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo $usage_str
            exit 1
            ;;
    esac
done

# Check mandatory arguments
if [[ -z "$instances_dir" ]]; then
    echo "Error: --instances-dir is mandatory."
    echo $usage_str
    exit 1
fi

if [[ -z "$instances_list" ]]; then
    echo "Error: --instances-list is mandatory."
    echo $usage_str
    exit 1
fi

if [[ -z "$configs_list" ]]; then
    echo "Error: --configs-list is mandatory."
    echo $usage_str
    exit 1
fi

if [[ -z "$time_limit" ]]; then
    echo "Error: --time-limit is mandatory."
    echo $usage_str
    exit 1
fi

if [[ -z "$num_seeds" ]]; then
    echo "Error: --num-seeds is mandatory."
    echo $usage_str
    exit 1
fi

# Function to resolve absolute path
resolve_path() {
    realpath "$1" 2>/dev/null
}

# Resolve paths to absolute and check for failure
instances_dir=$(resolve_path "$instances_dir")
if [[ -z "$instances_dir" ]]; then
    echo "Error: Failed to resolve path for --instances-dir."
    exit 1
fi

instances_list=$(resolve_path "$instances_list")
if [[ -z "$instances_list" ]]; then
    echo "Error: Failed to resolve path for --instances-list."
    exit 1
fi

solutions_dir=$(resolve_path "$solutions_dir")
mkdir $solutions_dir
if [[ -z "$solutions_dir" ]]; then
    echo "Error: Failed to resolve path for --solutions-dir."
    exit 1
fi
if [[ ! -d "$solutions_dir" ]]; then
    echo "Error: --solutions-dir '$solutions_dir' is not a valid directory."
    exit 1
fi

configs_dir=$(resolve_path "$configs_dir")
if [[ -z "$configs_dir" ]]; then
    echo "Error: Failed to resolve path for --configs-dir."
    exit 1
fi

configs_list=$(resolve_path "$configs_list")
if [[ -z "$configs_list" ]]; then
    echo "Error: Failed to resolve path for --configs-list."
    exit 1
fi

# Validate directories
if [[ ! -d "$instances_dir" ]]; then
    echo "Error: instances-dir '$instances_dir' is not a valid directory."
    exit 1
fi

if [[ ! -d "$configs_dir" ]]; then
    echo "Error: configs-dir '$configs_dir' is not a valid directory."
    exit 1
fi

# Validate files
if [[ ! -f "$instances_list" ]]; then
    echo "Error: instances-list '$instances_list' is not a valid file."
    exit 1
fi

if [[ ! -f "$configs_list" ]]; then
    echo "Error: configs-list '$configs_list' is not a valid file."
    exit 1
fi

# Validate time-limit
if ! [[ "$time_limit" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: time-limit must be a floating point number."
    exit 1
fi

# Validate num-seeds
if ! [[ "$num_seeds" =~ ^[0-9]+$ ]] || (( num_seeds <= 0 )); then
    echo "Error: --num-seeds must be an integer greater than 0."
    exit 1
fi

# Run experiments
./script_run_experiments.sh --instances-dir $instances_dir --instances-list $instances_list --solutions-dir $solutions_dir --configs-dir $configs_dir --configs-list $configs_list --percentage-random-sample $percentage_random_sample --time-limit $time_limit --num-seeds $num_seeds

file_list_instances_used="$solutions_dir/sorted_files_list.txt"
echo "Note: the list of instances used in the experiments was stored in the file ${file_list_instances_used}"
echo "This file should be the --instances-list arg value when generating the results by running script_generate_results.sh ." 

# Generate results
./script_generate_results.sh --instances-dir $instances_dir --instances-list $file_list_instances_used --solutions-dir $solutions_dir --configs-dir $configs_dir --configs-list $configs_list --time-limit $time_limit --num-seeds $num_seeds

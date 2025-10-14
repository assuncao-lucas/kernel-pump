#!/bin/bash

# Default values
# Mandatory arguments (no defaults)
instances_dir=""
instances_list=""
solutions_dir=""
configs_dir=""
configs_list=""
time_limit=""
num_seeds=""

# optional argument
percentage_random_sample=100

files=()
seeds=()
configs_vec=()

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
        --solutions-dir)
            solutions_dir="$2"
            shift 2
            ;;
        --configs-dir)
            configs_dir="$2"
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
            echo "Usage: $0 --instances-dir PATH --instances-list FILE --solutions-dir PATH --configs-dir PATH --configs-list FILE --time-limit FLOAT --num-seeds INT>0 [--percentage-random-sample INTEGER 0-100]"
            exit 1
            ;;
    esac
done

# Check mandatory arguments
if [[ -z "$instances_dir" ]]; then
    echo "Error: --instances-dir is mandatory."
    exit 1
fi

if [[ -z "$instances_list" ]]; then
    echo "Error: --instances-list is mandatory."
    exit 1
fi

if [[ -z "$solutions_dir" ]]; then
    echo "Error: --solutions-dir is mandatory."
    exit 1
fi

if [[ -z "$configs_dir" ]]; then
    echo "Error: --configs-dir is mandatory."
    exit 1
fi

if [[ -z "$configs_list" ]]; then
    echo "Error: --configs-list is mandatory."
    exit 1
fi

if [[ -z "$time_limit" ]]; then
    echo "Error: --time-limit is mandatory."
    exit 1
fi

if [[ -z "$num_seeds" ]]; then
    echo "Error: --num-seeds is mandatory."
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

# Absolute path where results will be stored
absolute_solutions_dir=$(realpath "$solutions_dir")

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

# Generate seeds sequence
seeds=($(seq 1 $num_seeds))

# Print results
echo "instances-dir: $instances_dir"
echo "instances-list: $instances_list"
echo "solutions-dir: $solutions_dir"
echo "configs-dir: $configs_dir"
echo "configs-list: $configs_list"
echo "percentage-random-sample: $percentage_random_sample"
echo "time-limit: $time_limit"
echo "num-seeds: $num_seeds"

# Move to build directory
cd ./build

# Build executable
echo "Building executable..."
make -j12

# ----------------------------
# Read and process instances
# ----------------------------
mapfile -t all_files < <(grep -v '^[[:space:]]*$' "$instances_list")


# Validate percentage_random_sample is integer between 0 and 100
if ! [[ "$percentage_random_sample" =~ ^[0-9]+$ ]] || (( percentage_random_sample < 0 || percentage_random_sample > 100 )); then
    echo "Error: --percentage-random-sample must be an integer between 0 and 100."
    exit 1
fi

if [[ "$percentage_random_sample" -lt 100 ]]; then
    total_files=${#all_files[@]}
    num_files_to_select=$(( total_files * percentage_random_sample / 100 ))

    # Ensure at least 1 file is selected
    if (( num_files_to_select < 1 )); then
        num_files_to_select=1
    fi

    # Shuffle all_files, select first num_files_to_select
    mapfile -t shuffled_files < <(printf '%s\n' "${all_files[@]}" | shuf)
    files=("${shuffled_files[@]:0:num_files_to_select}")

    # Sort the selected files alphabetically
    readarray -t sorted_files < <(printf '%s\n' "${files[@]}" | sort)

    echo "Randomly selected ${#sorted_files[@]} files ($percentage_random_sample% of total)"
else
    # Not random -> all files
    files=("${all_files[@]}")
    readarray -t sorted_files < <(printf '%s\n' "${files[@]}" | sort)
    echo "All ${#sorted_files[@]} files from instances_list are considered"
fi

# Print sorted_files for verification
echo "Files considered (sorted_files):"
for f in "${sorted_files[@]}"; do
    echo "  $f"
done

# ----------------------------
# Save the list of sorted files to a new file
sorted_files_list="$absolute_solutions_dir/sorted_files_list.txt"
printf '%s\n' "${sorted_files[@]}" > "$sorted_files_list"
echo "List of sorted files saved at: $sorted_files_list"


# ----------------------------
# Read configs
# ----------------------------
mapfile -t configs_vec < <(grep -v '^[[:space:]]*$' "$configs_list")

# ----------------------------
# Run Kernel Pump for all combinations of instances, configs and seeds
# ----------------------------
echo "Running tests for all combinations of instances, configs and seeds..."
for curr_instance in "${sorted_files[@]}"; do
    instance_base="${curr_instance%.mps}"       # Remove .mps extension
    for config_name in "${configs_vec[@]}"; do
        config_base="${config_name%.cfg}"      # Remove .cfg extension
        echo "* inst: $curr_instance"
        curr_instance_path="$instances_dir/$curr_instance"
        echo "* inst path: $curr_instance_path"
        config_path="$configs_dir/$config_name"
        echo "* config: $config_name"
        echo "* config path: $config_path"
        for seed in "${seeds[@]}"; do
            echo "- seed: $seed"
            ./kp $curr_instance_path --config $config_path seed=$seed solutionFolder=$absolute_solutions_dir timeLimit=$time_limit
            # Print absolute path of solution file
            sol_file="${absolute_solutions_dir}/s_${config_base}_${instance_base}_${seed}.sol"

            absolute_curr_solution_path=$(realpath "$sol_file")
            echo "Execution file generated at: $absolute_curr_solution_path"
        done
    done
done

# move back directory
cd ..
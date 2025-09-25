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
            echo "Usage: $0 --instances-dir PATH --instances-list FILE --solutions-dir PATH --configs-dir PATH --configs-list FILE --time-limit FLOAT --num-seeds INT>0"
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
echo "time-limit: $time_limit"
echo "num-seeds: $num_seeds"

# ----------------------------
# Read and process instances
# ----------------------------
mapfile -t all_files < "$instances_list"

# Not random -> all files
files=("${all_files[@]}")
readarray -t sorted_files < <(printf '%s\n' "${files[@]}" | sort)
echo "All ${#sorted_files[@]} files from instances_list are considered"

# Print sorted_files for verification
echo "Files considered (sorted_files):"
for f in "${sorted_files[@]}"; do
    echo "  $f"
done

# ----------------------------
# Read configs
# ----------------------------
cont=0
while IFS= read -r LINE_CONF; do
    configs_vec[$cont]="$LINE_CONF"
    ((cont++))
done < "$configs_list"

# ----------------------------
# Generate summary CSV
# ----------------------------
summary_csv="${absolute_solutions_dir}/summary.csv"
echo "Creating summary CSV..."
echo "instance,config,seed,status,time_kernel,total_time,iterations,buckets,last_bucket,first_bucket,value,reopt_value,real_gap,proj_gap,num_frac" > "$summary_csv"

for curr_instance in "${sorted_files[@]}"; do
    instance_base="${curr_instance%.mps}"       # Remove .mps extension
    for config_name in "${configs_vec[@]}"; do
        config_base="${config_name%.cfg}"      # Remove .cfg extension
        for seed in "${seeds[@]}"; do
            sol_file="${absolute_solutions_dir}/s_${config_base}_${instance_base}_${seed}.sol"

            if [[ ! -f "$sol_file" ]]; then
                echo "Warning: Solution file $sol_file not found. Skipping."
                continue
            fi

            status=$(grep -m1 "^STATUS:" "$sol_file" | cut -d':' -f2- | xargs)
            time_kernel=$(grep -m1 "^time building kernel" "$sol_file" | cut -d':' -f2- | xargs)
            total_time=$(grep -m1 "^total time" "$sol_file" | cut -d':' -f2- | xargs)

            # Enforce time limit
            if [[ -n "$total_time" ]]; then
                total_time_val=$(echo "$total_time" | awk '{print $1}')
                if (( $(echo "$total_time_val > $time_limit" | bc -l) )); then
                    total_time="$time_limit"
                fi
            fi

            iterations=$(grep -m1 "^# iterations" "$sol_file" | cut -d':' -f2- | xargs)
            buckets=$(grep -m1 "^# buckets" "$sol_file" | cut -d':' -f2- | xargs)
            last_bucket=$(grep -m1 "^last bucket visited" "$sol_file" | cut -d':' -f2- | xargs)
            first_bucket=$(grep -m1 "^first bucket to iter pump" "$sol_file" | cut -d':' -f2- | xargs)
            value=$(grep -m1 "^value" "$sol_file" | cut -d':' -f2- | xargs)
            reopt_value=$(grep -m1 "^reopt value" "$sol_file" | cut -d':' -f2- | xargs)
            real_gap=$(grep -m1 "^real integrality gap" "$sol_file" | cut -d':' -f2- | xargs)
            proj_gap=$(grep -m1 "^projection integrality gap" "$sol_file" | cut -d':' -f2- | xargs)
            num_frac=$(grep -m1 "^num frac" "$sol_file" | cut -d':' -f2- | xargs)

            echo "${instance_base},${config_base},${seed},${status},${time_kernel},${total_time},${iterations},${buckets},${last_bucket},${first_bucket},${value},${reopt_value},${real_gap},${proj_gap},${num_frac}" >> "$summary_csv"
        done
    done
done

# Print absolute path of summary CSV
absolute_csv_path=$(realpath "$summary_csv")
echo "Summary CSV generated at: $absolute_csv_path"

# ----------------------------
# Final step: Build and run results executable
# ----------------------------
echo "Entering ./results to build and run results executable..."
cd ./results/build

# Build results executable
echo "Building results executable..."
make -j12

# Reconstruct arguments (absolute paths, excluding random sampling flags)
results_args=( 
    --instances-dir "$instances_dir"
    --instances-list "$instances_list"
    --configs-list "$configs_list"
    --solutions-dir "$absolute_solutions_dir"
    --time-limit "$time_limit"
    --num-seeds "$num_seeds"
)

echo "Running results executable with arguments: ${results_args[*]}"
./results "${results_args[@]}"

# Go back to root directory
cd ../..

# ----------------------------
# Final message
# ----------------------------

echo "All raw solution files, the summary CSV, and the skeleton of the results LaTeX tables were stored in: $absolute_solutions_dir"
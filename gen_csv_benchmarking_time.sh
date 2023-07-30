#!/bin/bash

problems=("NonlinearOCP1b" "QuadcopterAD" "QuadcopterFull" "HangingChain" "MultiRTAC" "Nagumo")

# Define path to .csv file
path_to_output="${HOME}/Thesis_project/benchmarking_time_GN/"

# Define path to profiling tool command
path_to_profiling_tool="${HOME}/kokkos-tools/profiling/simple-kernel-timer/kp_reader"

output_file="output.csv"

# Remove existing output file
cd "$path_to_output"
rm -f "$output_file"

# Create an associative array to store the data for each problem
declare -A problem_data

for p in "${problems[@]}"; do

    # Initialize an empty string to store the data for the current problem
    problem_data[$p]=""

    # Go to the built application directory
    cd "${HOME}/Thesis_project/benchmarking_time_GN/$p"

    # Loop through .dat files starting with "filename"
    for file in nicholas-Inspiron-5567-*.dat; do
        # Extract the desired part from each file and append it to the data for the current problem
        data=$($path_to_profiling_tool "$file" | grep " -> Percentage in Kokkos kernels:" | sed 's/   -> Percentage in Kokkos kernels:                                   //')
        problem_data[$p]+="$data,"
    done

done

# Loop through the problems and write the problem name and its corresponding data in a new row
for p in "${problems[@]}"; do
    # Write the problem name and its corresponding data in a new row
    echo "$p,${problem_data[$p]}" >> "${path_to_output}${output_file}"
done

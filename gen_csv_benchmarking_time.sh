#!/bin/bash

problems=("NonlinearOCP1b" "QuadcopterAD" "QuadcopterFull" "HangingChain" "MultiRTAC")

# Define path to .csv file
path_to_output="/home/nicholas/Thesis_project/benchmarking_time/"

# Define path to profiling tool command
path_to_profiling_tool="/home/nicholas/kokkos-tools/profiling/simple-kernel-timer/kp_reader"

output_file="output.csv"

# Remove existing output file
cd $path_to_output
rm -f "$output_file"

for p in "${problems[@]}"
do

    echo "$p" >> "${path_to_output}${output_file}"

    # Go to the built application directory
    cd "/home/nicholas/Thesis_project/benchmarking_time/$p"

    # Loop through .dat files starting with "filename"
    for file in nicholas-Inspiron-5567-*.dat; do
        # Extract the desired part from each file and append it to the output file
        $path_to_profiling_tool "$file" | grep " -> Percentage in Kokkos kernels:" | sed 's/ -> Percentage in Kokkos kernels: //' >> "${path_to_output}${output_file}"
    done

done



#!/bin/bash

# Open the build directory and export kp_kernel_timer.so to the path
cd "/home/nicholas/Thesis_project/build/applications"
export KOKKOS_TOOLS_LIBS="/home/nicholas/kokkos-tools/profiling/simple-kernel-timer/kp_kernel_timer.so"

problems=("NonlinearOCP1b" "QuadcopterAD" "QuadcopterFull" "HangingChain" "MultiRTAC" "Nagumo")
headers=("nonlinear_dynamics.hpp" "quadcopter.hpp" "quadcopter_full.hpp" "hanging_chain.hpp" "multi-RTAC.hpp" "nagumo_schlogl.hpp") 

prev_element="${problems[0]}"
prev_header="${headers[0]}"

h="${headers[0]}"

j=0

for p in "${problems[@]}"
do

    cpp_file="/home/nicholas/Thesis_project/applications/tests_benchmarking.cpp"
    hpp_file="/home/nicholas/Thesis_project/problems/$h"

    # Modify cpp file to run benchmarking on the next problem 
    sed -i 's/#include <'"$prev_header"'>/#include <'"$h"'>/g' $cpp_file
    sed -i 's/alpaqa::TypeErasedControlProblem<config_t>::make<'"$prev_element"'>/alpaqa::TypeErasedControlProblem<config_t>::make<'"$p"'>/g' $cpp_file

    # Define the command to execute
    command="./tests_benchmarking"

    # Loop N times and run the command
    for i in {1..10}
    do
        sed -i 's/unsigned long int n_seed = '"$i"'/unsigned long int n_seed = '"$((i+1))"'/g' $hpp_file
        echo "Running command: $command for $p"
        
        # Build modified application/test
        cd /home/nicholas/Thesis_project/
        cmake --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ -S/home/nicholas/Thesis_project -B/home/nicholas/Thesis_project/build -G Ninja
        cmake --build build -j 

        # Go to the built application directory
        cd "/home/nicholas/Thesis_project/build/applications"
        
        $command

    done

    src_dir="/home/nicholas/Thesis_project/build/applications/"
    dest_dir="/home/nicholas/Thesis_project/benchmarking_time/$p"
    mkdir -p "$dest_dir"

    find "$src_dir" -type f -name 'nicholas-Inspiron-5567*' -print0 | while read -d '' file
    do
    # Check if the file is a dat file
    if [[ $(file -b "$file") == *"dat"* ]]; then
        # Move the file to the destination directory
        mv "$file" "$dest_dir"
    fi
    done

    prev_element="$p"
    prev_header="${headers[j]}"

    j=$((j+1))

    echo "j is $j"

    h="${headers[j]}"

done



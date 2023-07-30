#!/bin/bash

# Open the build directory and export kp_kernel_timer.so to the path
cd "${HOME}/Thesis_project/build/applications"

problems=("NonlinearOCP1b" "QuadcopterAD" "QuadcopterFull" "HangingChain" "MultiRTAC" "Nagumo")
headers=("nonlinear_dynamics.hpp" "quadcopter.hpp" "quadcopter_full.hpp" "hanging_chain.hpp" "multi-RTAC.hpp" "nagumo_schlogl.hpp") 

prev_element="${problems[0]}"
prev_header="${headers[0]}"

h="${headers[0]}"

j=0

for p in "${problems[@]}"
do

    cpp_file="${HOME}/Thesis_project/applications/test_tau_GN.cpp"
    hpp_file="${HOME}/Thesis_project/problems/$h"

    # Modify cpp file to run benchmarking on the next problem 
    sed -i 's/#include <'"$prev_header"'>/#include <'"$h"'>/g' $cpp_file
    sed -i 's/alpaqa::TypeErasedControlProblem<config_t>::make<'"$prev_element"'>/alpaqa::TypeErasedControlProblem<config_t>::make<'"$p"'>/g' $cpp_file
    sed -i 's/    std::string problem_name = '""$prev_element""'/    std::string problem_name = '""$p""'/g' $cpp_file
    # Define the command to execute
    command="./test_tau_GN"

    # Loop N times and run the command
    for i in {1..10}
    do
        sed -i 's/unsigned long int n_seed = '"$i"'/unsigned long int n_seed = '"$((i+1))"'/g' $hpp_file
        sed -i 's/    index_t sim_index = '"$i"'/    index_t sim_index = '"$((i+1))"'/g' $cpp_file
        echo "Running command: $command for $p"
        
        # Build modified application/test
        cd ${HOME}/Thesis_project/
        cmake --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ -S${HOME}/Thesis_project -B${HOME}/Thesis_project/build -G Ninja
        cmake --build build -j 

        # Go to the built application directory
        cd "${HOME}/Thesis_project/build/applications"
        
        $command

    done

    sed -i 's/    index_t sim_index = 11;/    index_t sim_index = 1;/g' $cpp_file

    src_dir="${HOME}/Thesis_project/build/applications/"
    dest_dir="${HOME}/Thesis_project/investigating_tau_GN/$p"
    mkdir -p "$dest_dir"

    find "$src_dir" -type f -name 'results_' -print0 | while read -d '' file
    do
    mv "$file" "$dest_dir"
    done

    prev_element="$p"
    prev_header="${headers[j]}"

    j=$((j+1))

    echo "j is $j"

    h="${headers[j]}"

done



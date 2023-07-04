#!/bin/bash

n_processor=("4" "3" "2")

prev_element="${n_processor[0]}"

hpp_file="${HOME}/Thesis_project/problems/multi-RTAC.hpp"

j=0

for p in "${n_processor[@]}"
do

    cd "${HOME}/Thesis_project/build/applications"
    
    cpp_file="${HOME}/Thesis_project/applications/tests_ss_vs_ms.cpp"

    # Modify cpp file to run benchmarking on the next problem 
    sed -i 's/const int np = '"$prev_element"';/const int np = '"$p"';/g' $cpp_file

    # Define the command to execute
    command="./tests_ss_vs_ms"

    # Loop N times and run the command
    for i in {2..4}
    do
        cd ${HOME}/Thesis_project/problems
        sed -i 's/"p_Nc = '"$i"',/p_Nc = '"$((i+1))"',/g' $hpp_file
        echo "Running command: $command for $p"
        
        # Build modified application/test
        cd ${HOME}/Thesis_project/
        cmake --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ -S${HOME}/Thesis_project -B${HOME}/Thesis_project/build -G Ninja
        cmake --build build -j 

        # Go to the built application directory
        cd "${HOME}/Thesis_project/build/applications"
        
        $command

    done

    prev_element="$p"

    j=$((j+1))

done



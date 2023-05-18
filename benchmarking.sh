#!/bin/bash

# Open the build directory and export kp_kernel_timer.so to the path
cd "/home/nicholas/Thesis_project/build/applications"
export KOKKOS_TOOLS_LIBS="/home/nicholas/kokkos-tools/profiling/simple-kernel-timer/kp_kernel_timer.so"

problems=("LinearOCP" 
        "Quadcopter" 
        "NonlinearOCP1")

prev_element="${problems[2]}"

for p in "${problems[@]}"
do

    cpp_file="/home/nicholas/Thesis_project/applications/tests_benchmarking.cpp"

    # Modify cpp file to run benchmarking on the next problem 
    sed -i 's/alpaqa::TypeErasedControlProblem<config_t>::make<'"$prev_element"'>/alpaqa::TypeErasedControlProblem<config_t>::make<'"$p"'>/g' $cpp_file

    # Build modified application/test
    cd /home/nicholas/Thesis_project/
    cmake --no-warn-unused-cli -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_BUILD_TYPE:STRING=Debug -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++ -S/home/nicholas/Thesis_project -B/home/nicholas/Thesis_project/build -G Ninja
    cmake --build build -j 

    # Go to the built application directory
    cd "/home/nicholas/Thesis_project/build/applications"

    # Define the command to execute
    command="./tests_benchmarking"

    # Loop N times and run the command
    for i in {1..3}
    do
        echo "Running command: $command for $p"
        $command
    done

    src_dir="/home/nicholas/Thesis_project/build/applications/"
    dest_dir="/home/nicholas/Thesis_project/build/applications/benchmarking/$p"
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

done



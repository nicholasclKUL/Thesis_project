find_package(casadi REQUIRED)
find_package(alpaqa REQUIRED)
find_package(Threads REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Wrapper to load the generated problem
add_library(problem-loader SHARED "problem-loader.cpp")
target_link_libraries(problem-loader PUBLIC alpaqa::casadi-loader)
target_link_libraries(problem-loader PUBLIC alpaqa::casadi-ocp-loader)
include(GenerateExportHeader)
generate_export_header(problem-loader
    EXPORT_FILE_NAME export/problem-loader-export.h)
target_include_directories(problem-loader 
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/export>)
set_target_properties(problem-loader PROPERTIES
    C_VISIBILITY_PRESET "hidden"
    VISIBILITY_INLINES_HIDDEN true)
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_options(problem-loader PRIVATE "LINKER:--exclude-libs,ALL")
endif()

# Uses Python to generate a CasADi problem
function(casadi_function_codegen_python target pythonfile formulation horizon)
    set(problem_name "${target}_${formulation}_${horizon}")
    add_custom_command(OUTPUT "${problem_name}.c"
                    COMMAND Python3::Interpreter
                    ARGS "${CMAKE_CURRENT_SOURCE_DIR}/${pythonfile}" "${target}"
                         "${formulation}" "${horizon}"
                    BYPRODUCTS "${problem_name}.tsv"
                    MAIN_DEPENDENCY "${pythonfile}")
    add_library("${target}" MODULE "${problem_name}.c")
    set_target_properties("${target}" PROPERTIES PREFIX "" DEBUG_POSTFIX ""
                          OUTPUT_NAME ${problem_name})
    add_custom_command(TARGET "${target}" POST_BUILD
                       COMMAND ${CMAKE_COMMAND}
                       ARGS -E copy "${problem_name}.tsv"
                            "$<TARGET_FILE_DIR:${target}>")
endfunction()

# Python CasADi problem selection
set(PROBLEM_NAME_OPTIONS "quadcopter" "hermans_bicycle" "hanging_chain")
set(PROBLEM_NAME "hanging_chain" CACHE STRING "Problem to benchmark")
set_property(CACHE PROBLEM_NAME PROPERTY STRINGS ${PROBLEM_NAME_OPTIONS})
if(NOT PROBLEM_NAME IN_LIST PROBLEM_NAME_OPTIONS)
    message(FATAL_ERROR "PROBLEM_NAME must be one of ${PROBLEM_NAME_OPTIONS}")
endif()
# Python CasADi problem formulation selection
set(PROBLEM_FORM_OPTIONS "ocp" "ss" "ss2")
set(PROBLEM_FORM "ss" CACHE STRING "Problem to benchmark")
set_property(CACHE PROBLEM_FORM PROPERTY STRINGS ${PROBLEM_FORM_OPTIONS})
if(NOT PROBLEM_FORM IN_LIST PROBLEM_FORM_OPTIONS)
    message(FATAL_ERROR "PROBLEM_FORM must be one of ${PROBLEM_FORM_OPTIONS}")
endif()
# Python CasADi problem horizon lenght selection
set(PROBLEM_HORIZON 40 CACHE STRING "Horizon length")

# Python CasADi problem generation
casadi_function_codegen_python(${PROBLEM_NAME}
    "../python/util/generate_problem.py"
    ${PROBLEM_FORM} ${PROBLEM_HORIZON})
target_compile_definitions(problem-loader PRIVATE
    PROBLEM_DLL=\"$<TARGET_FILE_NAME:${PROBLEM_NAME}>\"
    $<$<STREQUAL:${PROBLEM_FORM},ocp>:"CONTROL_PROBLEM">
    PROBLEM_SECOND_ORDER=$<IF:$<STREQUAL:${PROBLEM_FORM},ss2>,1,0>
    PROBLEM_HORIZON=${PROBLEM_HORIZON})
add_dependencies(problem-loader ${PROBLEM_NAME})

# Example C++ driver
add_executable(driver driver.cpp)
target_link_libraries(driver PRIVATE problem-loader)

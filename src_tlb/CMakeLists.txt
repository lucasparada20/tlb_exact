cmake_minimum_required(VERSION 3.2)
project(Tlb)

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Hard code the path to the CPLEX directory
set(CPLEX_DIR "/some/path/to/Cplex")
message("Using CPLEX_DIR = ${CPLEX_DIR}")

# Include directories for CPLEX
set(CPLEX_INCLUDE_DIRS 
    ${CPLEX_DIR}/cplex/include
    ${CPLEX_DIR}/concert/include
)

# Link directories for CPLEX
set(CPLEX_LIBRARY_DIRS 
    ${CPLEX_DIR}/cplex/lib/x86-64_linux/static_pic
    ${CPLEX_DIR}/concert/lib/x86-64_linux/static_pic
)

# Flags for CPLEX
set(CPLEX_FLAGS ilocplex concert cplex m pthread dl)

message("CPLEX_INCLUDE_DIRS = ${CPLEX_INCLUDE_DIRS}")
message("CPLEX_LIBRARY_DIRS = ${CPLEX_LIBRARY_DIRS}")
message("CPLEX_FLAGS = ${CPLEX_FLAGS}")

# Create a shared interface target for CPLEX 
add_library(CPLEX INTERFACE)
target_include_directories(CPLEX INTERFACE ${CPLEX_INCLUDE_DIRS})
target_link_directories(CPLEX INTERFACE ${CPLEX_LIBRARY_DIRS})
target_link_libraries(CPLEX INTERFACE ${CPLEX_FLAGS})

# Check that all required Cplex paths and libraries exist
foreach(dir IN LISTS CPLEX_INCLUDE_DIRS CPLEX_LIBRARY_DIRS)
    if(NOT EXISTS ${dir})
        message(FATAL_ERROR "Required CPLEX path not found: ${dir}")
    endif()
endforeach()

# === LEMON ===
set(LEMON_DIR "/some/path/to/lemon-1.3.1")
message("Using LEMON_DIR = ${LEMON_DIR}")

add_library(LEMON INTERFACE)
target_include_directories(LEMON INTERFACE ${LEMON_DIR})

# === Sources ===
# main_tlb_test.cpp
# main_tlb_train.cpp

set(SOURCES
	main_tlb_train.cpp
	LoadTLB.cpp
	RecourseCalculation.cpp
	ScenarioGraph.cpp
	Parameters.cpp
	McfpSolvers.cpp
	Network.cpp
	ParetoOpt.cpp
	ExactTLB.cpp
	ExactTLBCallBacks.cpp
	DP.cpp
)

option(USE_OMP "Use OMP to compute the recourse in parallel" OFF)

# === Always add and link OpenMP ===
find_package(OpenMP REQUIRED)
add_library(OpenMPFlags INTERFACE)
target_compile_options(OpenMPFlags INTERFACE ${OpenMP_CXX_FLAGS})
target_link_libraries(OpenMPFlags INTERFACE ${OpenMP_CXX_LIB_NAMES})

# Only define macro if enabled
if (USE_OMP)
    target_compile_definitions(OpenMPFlags INTERFACE USE_OMP)
    message(STATUS "OpenMP ENABLED (with USE_OMP macro)")
else()
    message(STATUS "OpenMP linked (but USE_OMP macro not defined)")
endif()

# === Executable and links ===
add_executable(exec_tlb ${SOURCES})
target_link_libraries(exec_tlb PRIVATE LEMON CPLEX OpenMPFlags)




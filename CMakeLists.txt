cmake_minimum_required(VERSION 3.27)
project(turyn_golay)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(MPI_C_COMPILER "/Users/wilfried/homebrew/bin/mpicc")
set(MPI_CXX_COMPILER "/Users/wilfried/homebrew/bin/mpicxx")
find_package(MPI REQUIRED)
# Find OpenMP package
# Find OpenMPI package
find_package(OpenMP REQUIRED)
set(ENV{OMP_NUM_THREADS} 8)
add_executable(turyn_golay main.cpp
        utilities.cpp)
target_link_libraries(turyn_golay PRIVATE MPI::MPI_CXX)
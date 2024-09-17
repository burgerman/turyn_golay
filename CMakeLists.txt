cmake_minimum_required(VERSION 3.27)
project(turyn_golay)

set(CMAKE_CXX_STANDARD 17)
set(MPI_C_COMPILER "/Users/wilfried/homebrew/bin/mpicc")
set(MPI_CXX_COMPILER "/Users/wilfried/homebrew/bin/mpicxx")
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
add_executable(turyn_golay main.cpp
        utilities.cpp)
target_link_libraries(turyn_golay PRIVATE MPI::MPI_CXX)
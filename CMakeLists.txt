cmake_minimum_required(VERSION 3.27)
project(turyn_golay)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# Find OpenMP package
find_package(OpenMP REQUIRED)
add_executable(turyn_golay main.cpp
        utilities.cpp)
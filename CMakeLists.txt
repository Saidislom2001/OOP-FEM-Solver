cmake_minimum_required(VERSION 3.14)
project(FEMSolver CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3
find_package(Eigen3 3.3 REQUIRED CONFIG)

add_executable(fem_solver main.cpp)
target_link_libraries(fem_solver PRIVATE Eigen3::Eigen)

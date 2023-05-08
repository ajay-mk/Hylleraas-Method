cmake_minimum_required(VERSION 3.25)
project(Hylleraas)

set(CMAKE_CXX_STANDARD 17)

#Prerequisites
find_package(Eigen3 REQUIRED)

add_executable(Hylleraas main.cpp helper_functions.cpp helper_functions.h)
target_link_libraries(Hylleraas Eigen3::Eigen)

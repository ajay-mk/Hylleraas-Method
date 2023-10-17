cmake_minimum_required(VERSION 3.25)
project(Hylleraas)

set(CMAKE_CXX_STANDARD 17)

# Multiprecision can be turned off or on with the USE_MULTIPRECISION option
option(USE_MULTIPRECISION "Use multiprecision" OFF)

#Prerequisites
find_package(Eigen3 REQUIRED)

add_executable(Hylleraas main.cpp helper_functions.cpp helper_functions.h)
target_link_libraries(Hylleraas Eigen3::Eigen)

if (USE_MULTIPRECISION)
    add_definitions(-DUSE_MULTIPRECISION)
    find_package(Boost REQUIRED)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries(Hylleraas LINK_PUBLIC ${Boost_LIBRARIES})
endif()
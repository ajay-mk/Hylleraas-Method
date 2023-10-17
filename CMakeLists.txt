cmake_minimum_required(VERSION 3.25)
project(Hylleraas)

set(CMAKE_CXX_STANDARD 17)

# Multiprecision can be turned off or on with the USE_MULTIPRECISION option
option(ENABLE_MULTIPRECISION "Use multiprecision" ON)
option(ENABLE_OPENMP "Use OpenMP" ON)

add_executable(Hylleraas main.cpp helper_functions.cpp helper_functions.h)

# Eigen3 is required
find_package(Eigen3 REQUIRED)
target_link_libraries(Hylleraas Eigen3::Eigen)

# OpenMP is optional

if (ENABLE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        add_definitions(-DUSE_OPENMP)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
        target_link_libraries(Hylleraas OpenMP::OpenMP_CXX)
        target_compile_definitions(Hylleraas PUBLIC USE_OPENMP)
    endif ()
endif ()


# Multiprecision is optional
if (ENABLE_MULTIPRECISION)
    add_definitions(-DUSE_MULTIPRECISION)
    find_package(Boost REQUIRED)
    include_directories(${Boost_INCLUDE_DIR})
    target_link_libraries(Hylleraas LINK_PUBLIC ${Boost_LIBRARIES})
endif ()
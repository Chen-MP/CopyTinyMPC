cmake_minimum_required(VERSION 3.16)
project(wheel_legged_demo)

set(CMAKE_CXX_FLAGS "-std=c++11")

set(CMAKE_BUILD_TYPE "Debug")

set(EXEC_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# --- Eigen3  ---
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/Eigen/Eigen)
# find_package(Eigen3 REQUIRED)
# set(EIGEN3_INCLUDE_DIR /usr/local/include/eigen3)
# include_directories(${EIGEN3_INCLUDE_DIR}/eigen3)
# --- matplotlibcpp.h  ---
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include/matplotlibcpp)
find_package(Python3 COMPONENTS Development)
include_directories(${Python3_INCLUDE_DIRS})

# add_definitions(-DEIGEN_STACK_ALLOCATION_LIMIT=0)

# ----------------------- test_tinympc ----------------
add_executable(tinympc  test/test_tinympc.cpp
                        common/types_constants.cpp
                        tinympc/discretization.cpp
                        tinympc/lqr_dp.cpp
                        tinympc/admm.cpp)
target_compile_definitions(tinympc PRIVATE USE_PENDULUM_CONFIG) 
target_link_libraries(tinympc ${Python3_LIBRARIES})

# ----------------------- test_cartpole ----------------
add_executable(test_cartpole    test/test_cartpole.cpp
                                common/types_constants.cpp
                                tinympc/discretization.cpp
                                tinympc/lqr_dp.cpp
                                tinympc/admm.cpp)
target_compile_definitions(test_cartpole PRIVATE USE_PENDULUM_CONFIG) 
target_link_libraries(test_cartpole ${Python3_LIBRARIES})
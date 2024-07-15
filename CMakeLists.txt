# This is the CMakeLists of the project 
# 
# 
cmake_minimum_required(VERSION 3.15)
project(DynaMapp VERSION 1.0 LANGUAGES CXX)

set(BUILD_TEST TRUE) 
# Set the build type if not already specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "release")
endif()

# Activate C++ 14 
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Eigen3 REQUIRED)
find_package(GTest  REQUIRED)
find_package(Boost  REQUIRED)
find_package(OpenMP REQUIRED)

if(BUILD_TEST)
    add_executable(test_types test/test_types.cpp)
    target_link_libraries(test_types GTest::GTest GTest::Main Eigen3::Eigen ${Boost_LIBRARIES})
endif()


#add_executable(kalman_filter DynaMapp/kalman_filter.cpp)
#target_link_libraries(kalman_filter GTest::GTest GTest::Main Eigen3::Eigen ${Boost_LIBRARIES})
 







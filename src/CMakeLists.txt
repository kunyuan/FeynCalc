cmake_minimum_required (VERSION 2.6)
project (FeynCalc)
#set(CMAKE_CXX_COMPILER icpc)
#set(CMAKE_CXX_COMPILER g++)
#set(CMAKE_CXX_COMPILER clang)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wuninitialized")
include_directories(${FeynCalc_SOURCE_DIR})
#message("source dir:" ${FeynmanSimulator_SOURCE_DIR})

file(GLOB_RECURSE SRCS *.cpp)
file(GLOB_RECURSE HDRS *.h)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/..") 
#crucial to output the .exe file to the root of project, and CLion compatible

ADD_EXECUTABLE(feyncalc.exe  ${SRCS} ${HDRS})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

#install (TARGETS simulator.exe DESTINATION ${PROJECT_SOURCE_DIR}/..)

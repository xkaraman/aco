cmake_minimum_required(VERSION 3.1)
project(AntSystem CXX)

set (CMAKE_CXX_STANDARD 11)
include_directories(include)
set(SRC_FILES
    src/main.cpp
    src/ACO.cpp
    src/AntSystemSimple.cpp
    src/Ant.cpp
    )
add_executable(antsystem ${SRC_FILES})

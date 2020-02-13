cmake_minimum_required(VERSION 2.8.3)
project(AntSystem CXX)

include_directories(include)
set(SRC_FILES
    src/main.cpp
    src/ACO.cpp
    )
add_executable(antsystem ${SRC_FILES})

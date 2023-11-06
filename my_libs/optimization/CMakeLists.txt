cmake_minimum_required(VERSION 3.25)
project(Optimization LANGUAGES CXX VERSION 0.1)
set(CMAKE_CXX_STANDARD 20)

add_library(
    optimization 
    src/jacobian.cpp 
    src/root.cpp
)

target_include_directories(
    optimization
    PUBLIC
    include
)
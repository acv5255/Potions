cmake_minimum_required(VERSION 3.25)
project( Ode LANGUAGES CXX VERSION 0.1)
set(CMAKE_CXX_STANDARD 20)

add_library(
    ode
    STATIC
    src/explicit.cpp
    src/implicit.cpp
    src/solve_ode.cpp
)

target_include_directories(
    ode
    PUBLIC
    include
    ../optimization/include
)
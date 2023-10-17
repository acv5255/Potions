cmake_minimum_required(VERSION 3.24)
set(LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

project(Potions VERSION 0.1)

add_subdirectory(external/armadillo-code)
add_subdirectory(external/nlopt)
add_subdirectory(external/yaml-cpp)
add_subdirectory(external/matplotlib-cpp)
add_subdirectory(external/Catch2)

add_subdirectory(src)
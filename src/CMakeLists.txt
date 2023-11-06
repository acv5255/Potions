cmake_minimum_required(VERSION 3.25)

add_executable(
    potions
    main.cpp
)

include_directories(
    ../external/armadillo-code/include/
    ../external/matplotlib-cpp/
    ../my_libs/optimization/include/
    ../external/yaml-cpp/include/
    ../include
)
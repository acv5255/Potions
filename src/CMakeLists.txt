cmake_minimum_required(VERSION 3.25)

add_executable(
    potions
    main.cpp
    equilibrium.cpp
    kinetic.cpp
    model_inputs.cpp
    run_options.cpp
    utils.cpp
)

include_directories(
    ../external/armadillo-code/include/
    ../external/matplotlib-cpp/
    ../my_libs/optimization/include/
    ../my_libs/ode/include/
    ../external/yaml-cpp/include/
    ../include
)
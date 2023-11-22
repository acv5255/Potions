cmake_minimum_required(VERSION 3.25)

add_executable(
    potions
    main.cpp
    solution.cpp
    equilibrium.cpp
    kinetic.cpp
    model_inputs.cpp
    run_options.cpp
    utils.cpp
    chemfile.cpp
    database.cpp
)

target_include_directories(
    potions
    PRIVATE
    ../external/matplotplusplus/source/
    ../build/external/matplotplusplus/source/matplot/
    ../external/armadillo-code/include/
    ../my_libs/optimization/include/
    ../my_libs/ode/include/
    ../external/yaml-cpp/include/
    ../include
)

target_link_libraries(
    potions
    PRIVATE
    armadillo
    yaml-cpp
    matplot
    nodesoup
    optimization
    ode
)

target_link_directories(
    potions
    PRIVATE
    ../build/external/yaml-cpp/
    ../build/external/matplotplusplus/source/matplot
    ../build/my_libs/optimization/
    ../build/my_libs/ode/
)

add_compile_options(
    # -g
    -O3
)

set_target_properties(
    potions PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ../../
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ../../
)
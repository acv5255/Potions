cmake_minimum_required(VERSION 3.25)

add_executable(
    potions
    main.cpp
    equilibrium.cpp
    kinetic.cpp
    model_inputs.cpp
    run_options.cpp
    utils.cpp
    chemfile.cpp
    database.cpp
    solution.cpp
)

include_directories(
    ../external/armadillo-code/include/
    ../external/matplotlib-cpp/
    ../my_libs/optimization/include/
    ../my_libs/ode/include/
    ../external/yaml-cpp/include/
    ../include
)

target_link_libraries(
    potions
    PUBLIC
    armadillo
    yaml-cpp
    optimization
    ode
)

target_link_directories(
    potions
    PUBLIC
    ../build/external/yaml-cpp/
    ../build/my_libs/optimization/
    ../build/my_libs/ode/
)

add_compile_options(
    -g
)

set_target_properties(
    potions PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ../../
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ../../
)
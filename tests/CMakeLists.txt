cmake_minimum_required(VERSION 3.25)

set(LANGUAGE CXX)
set(CMAKE_CXX_STANDARD 20)


add_executable(
    test_input 
    test_inputs.cpp
    ../src/model_inputs.cpp
    ../src/chemfile.cpp
    ../src/database.cpp
    ../src/utils.cpp
)

add_executable(
    test_optimization
    test_optimization.cpp
)

target_link_libraries(
    test_input
    PUBLIC
    Catch2::Catch2WithMain
    yaml-cpp
    armadillo
)

target_link_libraries(
    test_optimization
    PUBLIC
    Catch2::Catch2WithMain
    optimization
    armadillo
)

add_compile_options(
    -g
)

include_directories(../include)
include_directories(../external/yaml-cpp/include)

enable_testing()

include(CTest)
include(Catch)
catch_discover_tests(test_input)
catch_discover_tests(test_optimization)
cmake_minimum_required(VERSION 3.24)
set(LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

project(Potions VERSION 0.1)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

add_subdirectory(external/armadillo-code)
add_subdirectory(external/yaml-cpp)
add_subdirectory(external/matplotplusplus)
add_subdirectory(external/Catch2)

include_directories(
    external/matplotplusplus/source/
    build/external/matplotplusplus/source/matplot/
)

link_libraries(
    matplot
)

# My libraries
add_subdirectory(my_libs/optimization)
add_subdirectory(my_libs/ode)
add_subdirectory(tests)

enable_testing()

add_subdirectory(src)
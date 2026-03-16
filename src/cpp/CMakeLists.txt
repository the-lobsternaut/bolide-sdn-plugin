cmake_minimum_required(VERSION 3.14)
project(bolide_sdn_plugin CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header-only library (types.h + parser.h)
add_library(bolide INTERFACE)
target_include_directories(bolide INTERFACE include)

# Test binary
add_executable(test_bolide tests/test_bolide.cpp)
target_link_libraries(test_bolide bolide)

enable_testing()
add_test(NAME bolide COMMAND test_bolide)

add_executable(test_sixdof tests/test_sixdof.cpp)
target_include_directories(test_sixdof PRIVATE include)
add_test(NAME bolide_sixdof COMMAND test_sixdof)

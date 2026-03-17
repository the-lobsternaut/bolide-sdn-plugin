cmake_minimum_required(VERSION 3.14)
project(bolide_sdn_plugin CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Static library with physics/analysis implementations
add_library(bolide_lib STATIC src/bolide.cpp)
target_include_directories(bolide_lib PUBLIC include)

# Test binary
add_executable(test_bolide tests/test_bolide.cpp)
target_link_libraries(test_bolide bolide_lib)

enable_testing()
add_test(NAME bolide COMMAND test_bolide)

add_executable(test_sixdof tests/test_sixdof.cpp)
target_include_directories(test_sixdof PRIVATE include)
add_test(NAME bolide_sixdof COMMAND test_sixdof)

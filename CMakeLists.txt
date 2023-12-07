cmake_minimum_required(VERSION 3.8)
project(mc_traffic_optimization)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)

# Set the include directory
include_directories(${PROJECT_SOURCE_DIR}/src)

# Add subdirectories for the graph and qubo components
add_subdirectory(src/graph)
add_subdirectory(src/qubo)
add_subdirectory(src/anneal)
add_subdirectory(src/schwefel)

# Add the main executable
add_executable(mcto src/main.cpp)

# Link the graph and qubo libraries with the main executable
target_link_libraries(mcto graph qubo)

# Enable testing and add the tests subdirectory
enable_testing()
add_subdirectory(tests)
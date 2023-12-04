# tests/CMakeLists.txt

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Add test cpp files
add_executable(runUnitTests
    test_main.cpp
    graph_test.cpp
    qubo_test.cpp
)

# Link test executable against gtest & gtest_main and the parts of your project being tested
target_link_libraries(runUnitTests ${GTEST_LIBRARIES} pthread graph qubo)

# Add the tests to CMake's testing system
add_test(NAME runUnitTests COMMAND runUnitTests)
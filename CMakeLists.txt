# Find and include Google Test
find_package(GTest REQUIRED)

if(GTEST_FOUND)
    message(STATUS "Google Test found. Including Google Test in the build.")
    message(STATUS "Google Test include directories: ${GTEST_INCLUDE_DIRS}")
    message(STATUS "Google Test library: ${GTEST_LIBRARIES}")
    include_directories(${GTEST_INCLUDE_DIRS})
else()
    message(STATUS "Google Test not found. Tests will not be built.")
endif()

# Add your test sources here
add_executable(test_suite test/test1.cpp test/test2.cpp)

# Link the test executable with Google Test
target_link_libraries(test_suite GTest::GTest GTest::Main)
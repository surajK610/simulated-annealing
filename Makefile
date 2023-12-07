CXX      = g++

GTEST_INCLUDE_DIR = /gpfs/runtime/opt/googletest/1.8.0/include
GTEST_LIB_DIR = /gpfs/runtime/opt/googletest/1.8.0/lib

CXXFLAGS = -Wall -pedantic -O3 -fopenmp -I src/anneal -I src/schwefel -I $(GTEST_INCLUDE_DIR)

SRC_DIR := src
BUILD_DIR := build

MAIN_EXECUTABLE = $(BUILD_DIR)/run_sa

TEST_DIR := tests
TEST_SRCS := $(wildcard $(TEST_DIR)/*.cpp)

TEST_OBJS := $(patsubst $(TEST_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(TEST_SRCS))
SRC_SRCS := $(wildcard $(SRC_DIR)/*.cpp)
SRC_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_SRCS))

TEST_EXECUTABLE = $(BUILD_DIR)/tests

GTEST_FLAGS = -lgtest -lgtest_main -pthread

$(shell mkdir -p $(BUILD_DIR))

# Build the main executable
$(MAIN_EXECUTABLE): $(SRC_OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SRC_OBJS) -o $(MAIN_EXECUTABLE)

# Build the test executable
$(TEST_EXECUTABLE): $(TEST_OBJS) $(SRC_OBJS) $(GTEST_HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(TEST_OBJS) $(SRC_OBJS) -o $(TEST_EXECUTABLE) $(GTEST_FLAGS)

# Rule to compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile test source files
$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY : test
test: $(TEST_EXECUTABLE)
	./$(TEST_EXECUTABLE)

.PHONY : clean
clean :
	rm -f $(MAIN_EXECUTABLE) $(TEST_EXECUTABLE) $(BUILD_DIR)/*.o

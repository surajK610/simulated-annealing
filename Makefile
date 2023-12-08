CXX = g++
NVCC = nvcc

CXXFLAGS = -Wall -O3 -fopenmp -I src -I src/anneal -I src/schwefel
NVCCFLAGS = -O3
LDFLAGS = -L$(GTEST_LIB_DIR) -lgtest -lgtest_main -pthread -lcudart -L /usr/local/cuda/lib64

SRC_DIR := src
BUILD_DIR := build

MAIN_EXECUTABLE = $(BUILD_DIR)/main

CPP_SRCS := $(shell find $(SRC_DIR) -name '*.cpp')
CU_SRCS := $(shell find $(SRC_DIR) -name '*.cu')

CPP_OBJS := $(CPP_SRCS:.cpp=.o)
CU_OBJS := $(CU_SRCS:.cu=.o)

CPP_OBJS := $(addprefix $(BUILD_DIR)/,$(CPP_OBJS))
CU_OBJS := $(addprefix $(BUILD_DIR)/,$(CU_OBJS))

OBJS := $(CPP_OBJS) $(CU_OBJS)

.PHONY: all clean

all: $(MAIN_EXECUTABLE)


# Compile C++ source files
$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA source files
$(BUILD_DIR)/%.o: %.cu
	mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(MAIN_EXECUTABLE): $(OBJS)
	$(CXX) $(OBJS)  $(CXXFLAGS) $(LDFLAGS) -o $@

clean:
	rm -rf $(BUILD_DIR)
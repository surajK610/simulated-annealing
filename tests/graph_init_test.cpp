// tests/graph_init_test.cpp

#include "gtest/gtest.h"
#include "../src/graph/graph.h"

class TrafficGraphInitTest : public ::testing::Test {
protected:
    TrafficGraphInitTest graph;

    void SetUp() override {
        // Setup code (if needed)
    }

    void TearDown() override {
        // Teardown code (if needed)
    }
};

TEST_F(TrafficGraphInitTest, InitializePoints_WithValidBounds) {
    unsigned int numPoints = 10;
    unsigned int xBound = 100;
    unsigned int yBound = 100;

    graph.initializePoints(numPoints, xBound, yBound);

    ASSERT_EQ(graph.points.size(), numPoints); // Check if the correct number of points are initialized
    for (const auto& point : graph.points) {
        ASSERT_TRUE(point.x >= 0 && point.x <= xBound);
        ASSERT_TRUE(point.y >= 0 && point.y <= yBound);
    }
}

TEST_F(TrafficGraphInitTest, InitializePoints_WithZeroBounds) {
    graph.initializePoints(10, 0, 0);
    ASSERT_TRUE(graph.points.empty()); // Expect no points to be initialized
}

TEST_F(TrafficGraphInitTest, InitializeGraph_WithValidParameters) {
    unsigned int numPoints = 10;
    unsigned int additionalEdges = 5;
    unsigned int xBound = 100;
    unsigned int yBound = 100;

    graph.initializeGraph(numPoints, additionalEdges, xBound, yBound);

    // Add checks for valid graph initialization, e.g., number of edges, connectivity, etc.
}

TEST_F(TrafficGraphInitTest, InitializeGraph_WithInsufficientPoints) {
    graph.initializeGraph(1, 5, 100, 100); // Only one point
    ASSERT_TRUE(graph.edges.empty()); // Expect no edges to be formed
}

// Add more tests here to cover different scenarios and edge cases


// Additional tests can include:
// - InitializeGraph with zero points
// - InitializeGraph with a large number of points
// - InitializePoints with large bounds
// - InitializePoints with small bounds
// - InitializeGraph with a large number of additional edges
// - InitializeGraph with zero additional edges

// int main(int argc, char **argv) {
//     ::testing::InitGoogleTest(&argc, argv);
//     return RUN_ALL_TESTS();
// }

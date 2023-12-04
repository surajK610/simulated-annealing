// tests/graph_test.cpp

#include "gtest/gtest.h"
#include "../src/graph/graph.h" // Update the path to your trafficGraph.h

// Test Fixture for setting up a TrafficGraph with Points and Edges
class TrafficGraphTest : public ::testing::Test {
protected:
    TrafficGraph graph;
    std::vector<Point> points;
    std::vector<Edge> edges;

    void SetUp() override {
        // Initialize Points
        points.push_back({0, 0});
        points.push_back({1, 1});
        // ... Add more points as needed

        // Initialize Edges
        Edge edge1 = {&points[0], &points[1], 1.0};
        // ... Add more edges as needed

        // Add points and edges to graph
        graph.points = points;
        graph.edges.push_back(edge1);
        // ... Add more edges to graph as needed
    }
};

// Test case for mapping endpoints to a car
TEST_F(TrafficGraphTest, MapEndpoints) {
    Car car;
    car.id = 1;
    // Assuming mapEndpoints assigns the first and last points as source and destination
    graph.mapEndpoints(car);

    EXPECT_EQ(car.source, &points.front());
    EXPECT_EQ(car.destination, &points.back());
}

// Test case for finding all paths
TEST_F(TrafficGraphTest, FindAllPaths) {
    auto paths = graph.findAllPaths(&points.front(), &points.back());
    // Assuming findAllPaths should find at least one path
    EXPECT_GT(paths.size(), 0);
}

// Test case for finding alternative paths
TEST_F(TrafficGraphTest, FindAlternativePaths) {
    std::vector<Route> allPaths = graph.findAllPaths(&points.front(), &points.back());
    auto alternativePaths = graph.findAlternativePaths(allPaths);
    
    // Assuming findAlternativePaths should return a subset of allPaths
    EXPECT_LE(alternativePaths.size(), allPaths.size());
}

// Add more tests as needed...

// The main function to run all tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
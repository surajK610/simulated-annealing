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
        points.push_back({2, 2}); // Added an extra point for the shortest path test
        // ... Add more points as needed

        // Initialize Edges
        Edge edge1 = {&points[0], &points[1], 1.0};
        Edge edge2 = {&points[1], &points[2], 2.0};
        Edge edge3 = {&points[0], &points[2], 3.0}; // Direct but longer path
        // ... Add more edges as needed
        edges.push_back(edge1);
        edges.push_back(edge2);
        edges.push_back(edge3);
        // Add points and edges to graph
        graph.points = points;
        graph.edges = edges;
        // ... Add more edges to graph as needed
    }
};

// Test case for Dijkstra's algorithm
TEST_F(TrafficGraphTest, FindShortestPath) {
    Route* shortestRoute = graph.findShortestPath(&points[0], &points[2]);

    ASSERT_NE(shortestRoute, nullptr);
    ASSERT_EQ(shortestRoute->pathLen, 2);
    EXPECT_EQ(shortestRoute->route[0].start, &points[0]);
    EXPECT_EQ(shortestRoute->route[0].end, &points[1]);
    EXPECT_EQ(shortestRoute->route[1].start, &points[1]);
    EXPECT_EQ(shortestRoute->route[1].end, &points[2]);

    delete shortestRoute;
}


TEST_F(TrafficGraphTest, FindAllPaths) {
    auto paths = graph.findAllPaths(&points.front(), &points.back());
    // Assuming findAllPaths should find at least one path
    EXPECT_GT(paths.size(), 0);
}

// Add more tests as needed...

// The main function to run all tests
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
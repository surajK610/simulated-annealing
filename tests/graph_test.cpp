// tests/graph_test.cpp

#include "gtest/gtest.h"
#include "../src/graph/graph.h"

class TrafficGraphTest : public ::testing::Test {
protected:
    TrafficGraph graph;
    std::vector<Point> points;
    std::vector<Edge> edges;

public:
    TrafficGraphTest() {
    }

    void SetUp() override {
        std::vector<Point> points;
        std::vector<Edge> edges;
        points.push_back({0, 0});
        points.push_back({1, 1});
        points.push_back({2, 2}); 

        Edge edge1 = {&points[0], &points[1], 1.0};
        Edge edge2 = {&points[1], &points[2], 2.0};
        Edge edge3 = {&points[0], &points[2], 3.0};

        edges.push_back(edge1);
        edges.push_back(edge2);
        edges.push_back(edge3);
        graph = TrafficGraph(points, edges);
    }

};

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

TEST_F(TrafficGraphTest, InitializePoints) {
    unsigned int numPoints = 10;
    unsigned int xBound = 100;
    unsigned int yBound = 100;

    graph.initializePoints(numPoints, xBound, yBound);

    ASSERT_EQ(graph.getPoints().size(), numPoints); // Check if the correct number of points are initialized
    for (const auto& point : graph.getPoints()) {
        ASSERT_TRUE(point.x >= 0 && point.x <= xBound);
        ASSERT_TRUE(point.y >= 0 && point.y <= yBound);
    }

    graph.initializePoints(10, 0, 0);
    ASSERT_TRUE(graph.getPoints().empty()); 

    graph.initializeGraph(1, 5, 100, 100); // Only one point
    ASSERT_TRUE(graph.getEdges().empty()); // Expect no edges to be formed
}

TEST_F(TrafficGraphTest, FindAllPaths) {
    auto paths = graph.findAllPaths(&points.front(), &points.back());
    EXPECT_GT(paths.size(), 0);
}
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
        points.push_back({0, 0});
        points.push_back({1, 1});
        points.push_back({2, 2});
        points.push_back({3, 3}); 


        Edge edge1 = {&points[0], &points[1], 1.0};
        Edge edge2 = {&points[1], &points[2], 2.0};
        Edge edge3 = {&points[2], &points[3], 2.0};

        edges.push_back(edge1);
        edges.push_back(edge2);
        edges.push_back(edge3);

        graph = TrafficGraph(&points, &edges);

    }

    void printPoints(const std::vector<Point>& points) {
        std::cout << "Points:" << std::endl;
        for (const auto& point : points) {
            std::cout << &point <<  " (" << point.x << ", " << point.y << ")" << std::endl;
        }
    }

    void printEdges(const std::vector<Edge>& edges) {
        std::cout << "Edges:" << std::endl;
        for (const auto& edge : edges) {
            std::cout << "Start: " << edge.start << " (" << edge.start->x << ", " << edge.start->y << "), "
                      << "End: " << edge.end << " (" <<  edge.end->x << ", " << edge.end->y << "), "
                      << "Distance: " << edge.distance << std::endl;
        }
    }

    void printRoute(const Route& route) {
        std::cerr << "Route length: " << route.pathLen << "\n";
        for (int i = 0; i < route.pathLen; ++i) {
            const auto& edge = route.route[i];
            std::cerr << "Edge " << i << ": (" << edge.start->x << ", " << edge.start->y << ") to ("
                      << edge.end->x << ", " << edge.end->y << "), Dist: " << edge.distance << "\n";
        }
    }

    void printAllPaths(const std::vector<Route>& paths) {
        std::cerr << "Total paths found: " << paths.size() << "\n";
        for (const auto& path : paths) {
            printRoute(path);
        }
    }

};


TEST_F(TrafficGraphTest, InitializePoints) {
    unsigned int numPoints = 10;
    unsigned int xBound = 100;
    unsigned int yBound = 100;

    graph.initializePoints(numPoints, xBound, yBound);

    std::cout << "Checking number of points initialized..." << std::endl;
    ASSERT_EQ(graph.getPoints().size(), numPoints);
    // printPoints(graph.getPoints());    
    
    for (const auto& point : graph.getPoints()) {
        ASSERT_TRUE(point.x >= 0 && point.x <= xBound);
        ASSERT_TRUE(point.y >= 0 && point.y <= yBound);
    }
    // std::err << "Test Case 0: Zero Bound" << std::endl;
    // graph.initializePoints(10, 0, 0);
    // ASSERT_TRUE(graph.getPoints().empty()); 
    // printPoints(graph.getPoints());    

    // graph.initializeGraph(1, 5, 100, 100); // Only one point
    // ASSERT_TRUE(graph.getEdges().empty()); // Expect no edges to be formed

    ////////////////////////////////////////////
    // Test Case 1: Initialize graph with a small number of points and edges
    std::cerr << "    ////////////////////////////////////////////" << std::endl;
    std::cerr << "    Test Case 1: Small number of points and edges" << std::endl;
    std::cerr << "    ////////////////////////////////////////////" << std::endl;

    graph.initializeGraph(5, 3, 50, 50);  // 5 points, 3 additional edges
    // printPoints(graph.getPoints());
    printEdges(graph.getEdges());

    ASSERT_FALSE(graph.getPoints().empty());
    ASSERT_FALSE(graph.getEdges().empty());

    // Test Case 2: Initialize graph with no additional edges
    std::cerr << "    ////////////////////////////////////////////" << std::endl;
    std::cerr << "    \nTest Case 2: No additional edges" << std::endl;
    std::cerr << "    ////////////////////////////////////////////" << std::endl;

    graph.initializeGraph(4, 0, 30, 30);  // 4 points, no additional edges
    // printPoints(graph.getPoints());
    // printEdges(graph.getEdges());

    ASSERT_FALSE(graph.getPoints().empty());
    //ASSERT_EQ(graph.getEdges().size(), 3);  // Edges forming a Minimum Spanning Tree

    // Test Case 3: Large graph
    std::cerr << "    ////////////////////////////////////////////" << std::endl;
    std::cerr << "Test Case 3: Large graph" << std::endl;
    std::cerr << "    ////////////////////////////////////////////" << std::endl;
    graph.initializeGraph(20, 10, 100, 100);  // 20 points, 10 additional edges
    // printPoints(graph.getPoints());
    // printEdges(graph.getEdges());

    ASSERT_FALSE(graph.getPoints().empty());
    ASSERT_FALSE(graph.getEdges().empty());

    // Test Case 4: Edge case with only 1 point (no edges should be formed)
    std::cerr << "Test Case 4: Single point" << std::endl;
    graph.initializeGraph(1, 5, 100, 100);  // 1 point, request for additional edges
    // printPoints(graph.getPoints());
    // printEdges(graph.getEdges());

    ASSERT_EQ(graph.getPoints().size(), 1);
    ASSERT_TRUE(graph.getEdges().empty());  // No edges possible with a single point

    // Test Case 5: Edge case with zero points
    std::cerr << "    ////////////////////////////////////////////" << std::endl;
    std::cerr << "Test Case 5: Zero points" << std::endl;
    std::cerr << "    ////////////////////////////////////////////" << std::endl;

    graph.initializeGraph(0, 5, 100, 100);  // 0 points, irrelevant number of additional edges
    printPoints(graph.getPoints());
    printEdges(graph.getEdges());
    std::cerr << "Edges: " << graph.getEdges().size() << std::endl;

    ASSERT_EQ(graph.getPoints().size(), 0);
    ASSERT_EQ(graph.getEdges().size(), 0);
}

TEST_F(TrafficGraphTest, FindShortestPath) {
    Route* shortestRoute = graph.findShortestPath(&points[0], &points[3]);

    std::cerr << "    ////////////////////////////////////////////" << std::endl;
    std::cerr << "Test Case 6: Shortest Path Unique" << std::endl;
    std::cerr << "    ////////////////////////////////////////////" << std::endl;

    ASSERT_NE(shortestRoute, nullptr);
    ASSERT_EQ(shortestRoute->pathLen, 3);
    EXPECT_EQ(shortestRoute->route[0].start, &points[0]);
    EXPECT_EQ(shortestRoute->route[0].end, &points[1]);
    EXPECT_EQ(shortestRoute->route[1].start, &points[1]);
    EXPECT_EQ(shortestRoute->route[1].end, &points[2]);
    EXPECT_EQ(shortestRoute->route[2].start, &points[2]);
    EXPECT_EQ(shortestRoute->route[2].end, &points[3]);
    delete shortestRoute;

    edges.clear();
    Edge edge1 = {&points[0], &points[1], 1.0};
    Edge edge2 = {&points[1], &points[2], 2.0};
    Edge edge3 = {&points[0], &points[2], 4.0};

    edges.push_back(edge1);
    edges.push_back(edge2);
    edges.push_back(edge3);

    graph = TrafficGraph(&points, &edges);
    shortestRoute = graph.findShortestPath(&points[0], &points[2]);

    std::cerr << "    ////////////////////////////////////////////" << std::endl;
    std::cerr << "Test Case 6: Shortest Path PathLen Greater" << std::endl;
    std::cerr << "    ////////////////////////////////////////////" << std::endl;

    ASSERT_NE(shortestRoute, nullptr);
    ASSERT_EQ(shortestRoute->pathLen, 2);
    EXPECT_EQ(shortestRoute->route[0].start, &points[0]);
    EXPECT_EQ(shortestRoute->route[0].end, &points[1]);
    EXPECT_EQ(shortestRoute->route[1].start, &points[1]);
    EXPECT_EQ(shortestRoute->route[1].end, &points[2]);
    delete shortestRoute;
}

TEST_F(TrafficGraphTest, FindAllPathsStart ) {
    std::cerr << "////////// Test: START //////////\n";
}


TEST_F(TrafficGraphTest, FindAllPaths_MediumGraph) {
    std::cerr << "////////// Test: FindAllPaths_MediumGraph //////////\n";
    // Define a medium complexity graph
    // Initialize the graph with predefined points and edges
    // Call findAllPaths and print the selected start and endpoint using cerr

    // Example
    graph.initializeGraph(5, 3, 100, 100);  // Five points, a few additional edges
    
    auto startPoint = &graph.getPoints().front();
    auto endPoint = &graph.getPoints().back();
    
    std::cerr << "Selected Start Point: (" << startPoint->x << "," << startPoint->y << ")\n";
    std::cerr << "Selected End Point: (" << endPoint->x << "," << endPoint->y << ")\n";

    auto paths = graph.findAllPaths(startPoint, endPoint);
    EXPECT_GT(paths.size(), 0);
    std::cerr << "Paths found: " << paths.size() << "\n";
    printAllPaths(paths);
}

TEST_F(TrafficGraphTest, FindAllPaths) {
    auto paths = graph.findAllPaths(&points.front(), &points.back());
    EXPECT_GT(paths.size(), 0);
}


TEST_F(TrafficGraphTest, FindAllPaths_ComplexGraph) {
    std::cerr << "////////// Test: FindAllPaths_ComplexGraph //////////\n";
    // Define a complex graph
    // Initialize the graph with predefined points and edges
    // Call findAllPaths and print the selected start and endpoint using cerr

    // Example
    std::cerr.setstate(std::ios_base::failbit);
    graph.initializeGraph(10, 5, 100, 100);  // Ten points, several additional edges
    std::cerr.clear();

    auto startPoint = &graph.getPoints().front();
    auto endPoint = &graph.getPoints().back();
    
    std::cerr << "Selected Start Point: (" << startPoint->x << "," << startPoint->y << ")\n";
    std::cerr << "Selected End Point: (" << endPoint->x << "," << endPoint->y << ")\n";

    auto paths = graph.findAllPaths(startPoint, endPoint);
    EXPECT_GT(paths.size(), 0);
    std::cerr << "Paths found: " << paths.size() << "\n";
    printAllPaths(paths);
}

TEST_F(TrafficGraphTest, FindAllPaths_RandomGraph) {
    std::cerr << "////////// Test: FindAllPaths_RandomGraph //////////\n";
    // Generate a random graph
    // Initialize the graph with random points and edges
    // Call findAllPaths and print the selected start and endpoint using cerr

    // Example
    std::cerr.setstate(std::ios_base::failbit);
    graph.initializeGraph(20, 10, 1000, 1000);  // Larger, randomly generated graph
    std::cerr.clear();

    auto startPoint = &graph.getPoints().front();
    auto endPoint = &graph.getPoints().back();
    
    std::cerr << "Selected Start Point: (" << startPoint->x << "," << startPoint->y << ")\n";
    std::cerr << "Selected End Point: (" << endPoint->x << "," << endPoint->y << ")\n";

    auto paths = graph.findAllPaths(startPoint, endPoint);
    EXPECT_GT(paths.size(), 0);
    std::cerr << "Paths found: " << paths.size() << "\n";
    printAllPaths(paths);
}
// TEST_F(TrafficGraphTest, FindAllPaths) {
//     auto paths = graph.findAllPaths(&points.front(), &points.back());
//     EXPECT_GT(paths.size(), 0);
// }

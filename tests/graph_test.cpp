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
        Edge edge4 = {&points[0], &points[3], 12.0};

        edges.push_back(edge1);
        edges.push_back(edge2);
        edges.push_back(edge3);
        edges.push_back(edge4);

        graph = TrafficGraph(&points, &edges);

    }

    void SetUpLinearGraph() {
        points.clear();
        edges.clear();

        points.push_back({0, 0});
        points.push_back({1, 0});
        points.push_back({2, 0});
        points.push_back({3, 0}); 

        Edge edgeAB = {&points[0], &points[1], 1.0};
        Edge edgeBC = {&points[1], &points[2], 1.0};
        Edge edgeCD = {&points[2], &points[3], 1.0};

        edges.push_back(edgeAB);
        edges.push_back(edgeBC);
        edges.push_back(edgeCD);

        graph = TrafficGraph(&points, &edges);
    }

    void SetUpLoopGraph() {
        points.clear();
        edges.clear();

        points.push_back({0, 0});
        points.push_back({1, 0});
        points.push_back({1, 1});
        points.push_back({0, 1});

        Edge edgeAB = {&points[0], &points[1], 1.0};
        Edge edgeBC = {&points[1], &points[2], 1.0};
        Edge edgeCD = {&points[2], &points[3], 1.0};
        Edge edgeDA = {&points[3], &points[0], 1.0};

        edges.push_back(edgeAB);
        edges.push_back(edgeBC);
        edges.push_back(edgeCD);
        edges.push_back(edgeDA);

        graph = TrafficGraph(&points, &edges);
    }

    void SetUpRepeatEdgeGraph() {
        points.clear();
        edges.clear();

        points.push_back({0, 0}); 
        points.push_back({1, 0}); 

        Edge edge1 = {&points[0], &points[1], 1.0};
        Edge edge2 = {&points[0], &points[1], 2.0};
        Edge edge3 = {&points[0], &points[1], 3.0};

        edges.push_back(edge1);
        edges.push_back(edge2);
        edges.push_back(edge3);

        graph = TrafficGraph(&points, &edges);
    }

    void SetUpIsolatedPointsGraph() {
        points.clear();
        edges.clear();

        points.push_back({0, 0}); 
        points.push_back({1, 0}); 
        points.push_back({2, 0});
        points.push_back({3, 0});

        Edge edgeAB = {&points[0], &points[1], 1.0};

        edges.push_back(edgeAB);

        graph = TrafficGraph(&points, &edges);
    }

    void printPoints(const std::vector<Point>* points) {
        std::cout << "Points:" << std::endl;
        for (const auto& point : *points) {
            std::cout << &point <<  " (" << point.x << ", " << point.y << ")" << std::endl;
        }
    }

    void printEdges(const std::vector<Edge>* edges) {
        std::cout << "Edges:" << std::endl;
        for (const auto& edge : *edges) {
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

    // std::cout << "Checking number of points initialized..." << std::endl;
    ASSERT_EQ(graph.getPoints()->size(), numPoints);
    // printPoints(graph.getPoints());    
    
    for (const auto& point : *graph.getPoints()) {
        ASSERT_TRUE(point.x >= 0 && point.x <= xBound);
        ASSERT_TRUE(point.y >= 0 && point.y <= yBound);
    }
   
    // Test Case 1: Initialize graph with a small number of points and edges
    // std::cerr << "////////// Test Case 1: Small number of points and edges //////////\n";

    graph.initializeGraph(5, 3, 50, 50);  // 5 points, 3 additional edges
    // printPoints(graph.getPoints());
    printEdges(graph.getEdges());

    ASSERT_FALSE(graph.getPoints()->empty());
    ASSERT_FALSE(graph.getEdges()->empty());

    // Test Case 2: Initialize graph with no additional edges
    // std::cerr << "////////// Test Case 2: No additional edges //////////\n";


    graph.initializeGraph(4, 0, 30, 30);  // 4 points, no additional edges


    ASSERT_FALSE(graph.getPoints()->empty());
    ASSERT_EQ(graph.getEdges()->size(), 3);  // Edges forming a Minimum Spanning Tree

    // Test Case 3: Large graph
    // std::cerr << "////////// Test Case 3: Large graph //////////\n";

    graph.initializeGraph(20, 10, 100, 100);  // 20 points, 10 additional edges

    ASSERT_FALSE(graph.getPoints()->empty());
    ASSERT_FALSE(graph.getEdges()->empty());

    // Test Case 4: Edge case with only 1 point (no edges should be formed)
    // std::cerr << "////////// Test Case 4: Single point //////////\n";

    graph.initializeGraph(1, 5, 100, 100);  // 1 point, request for additional edges

    ASSERT_EQ(graph.getPoints()->size(), 1);
    ASSERT_TRUE(graph.getEdges()->empty());  // No edges possible with a single point

    // Test Case 5: Edge case with zero points
    // std::cerr << "////////// Test Case 5: Zero Points //////////\n";

    graph.initializeGraph(0, 5, 100, 100);  // 0 points, irrelevant number of additional edges
    printPoints(graph.getPoints());
    printEdges(graph.getEdges());
    std::cerr << "Edges: " << graph.getEdges()->size() << std::endl;

    ASSERT_EQ(graph.getPoints()->size(), 0);
    ASSERT_EQ(graph.getEdges()->size(), 0);
}

TEST_F(TrafficGraphTest, FindShortestPath) {
    Route* shortestRoute = graph.findShortestPath(&points[0], &points[3]);

    // std::cerr << "////////// Test Case 6: Shortest Path Unique //////////\n";

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

    // std::cerr << "////////// Test Case 6: Shortest Path PathLen Greater //////////\n";

    ASSERT_NE(shortestRoute, nullptr);
    ASSERT_EQ(shortestRoute->pathLen, 2);
    EXPECT_EQ(shortestRoute->route[0].start, &points[0]);
    EXPECT_EQ(shortestRoute->route[0].end, &points[1]);
    EXPECT_EQ(shortestRoute->route[1].start, &points[1]);
    EXPECT_EQ(shortestRoute->route[1].end, &points[2]);
    delete shortestRoute;
}

TEST_F(TrafficGraphTest, FindAllPathsStart ) {
    // std::cerr << "////////// Test: FindAllPathsStart //////////\n";
    std::vector<Route> allPaths = graph.findAllPaths(&points[0], &points[3]);
    ASSERT_EQ(allPaths.size(), 2);
    allPaths = graph.findAllPaths(&points[0], &points[2]);
    ASSERT_EQ(allPaths.size(), 2);
    allPaths= graph.findAllPaths(&points.front(), &points.back());
    EXPECT_GT(allPaths.size(), 0);
}


TEST_F(TrafficGraphTest, FindAllPathsMediumGraph) {
    // std::cerr << "////////// Test: FindAllPathsMediumGraph //////////\n";
    graph.initializeGraph(5, 3, 100, 100);  // Five points, a few additional edges
    
    auto startPoint = &graph.getPoints()->front();
    auto endPoint = &graph.getPoints()->back();

    auto paths = graph.findAllPaths(startPoint, endPoint);
    EXPECT_GT(paths.size(), 0);
    printAllPaths(paths);
}

TEST_F(TrafficGraphTest, FindAllPathsLinearGraph) {
    // std::cerr << "////////// Test: FindAllPathsLinearGraph //////////\n";
    SetUpLinearGraph(); 

    std::vector<Route> allPathsABCD = graph.findAllPaths(&points[0], &points[3]);
    ASSERT_EQ(allPathsABCD.size(), 1);  

    std::vector<Route> allPathsBC = graph.findAllPaths(&points[1], &points[2]);
    ASSERT_EQ(allPathsBC.size(), 1);  
}

TEST_F(TrafficGraphTest, FindAllPathsLoopGraph) {
    // std::cerr << "////////// Test: FindAllPathsLoopGraph //////////\n";
    SetUpLoopGraph();
    std::vector<Route> allPathsABDA = graph.findAllPaths(&points[0], &points[3]);
    ASSERT_EQ(allPathsABDA.size(), 2); 
}

TEST_F(TrafficGraphTest, FindAllPathsIsolatedPoints) {
    // std::cerr << "////////// Test: FindAllPathsIsolatedPoints //////////\n";
    SetUpIsolatedPointsGraph();

    std::vector<Route> allPathsCD = graph.findAllPaths(&points[2], &points[3]);
    ASSERT_EQ(allPathsCD.size(), 0); 

    std::vector<Route> allPathsAD = graph.findAllPaths(&points[0], &points[3]);
    ASSERT_EQ(allPathsAD.size(), 0);
}

TEST_F(TrafficGraphTest, FindAlternativePathsWithNoAlternatives) {
    SetUpLinearGraph();

    Point* source = &points[0];
    Point* destination = &points.back();

    std::vector<Route> alternativeRoutes = graph.findAlternativePaths(source, destination);
    ASSERT_TRUE(alternativeRoutes.empty());
}

TEST_F(TrafficGraphTest, FindAlternativePathsLoop) {
    SetUpLoopGraph();
    Point* source = &points[0];
    Point* destination = &points.back();
    std::vector<Route> alternativeRoutes = graph.findAlternativePaths(source, destination);
    std::cerr << "Alternative routes: " << alternativeRoutes.size() << std::endl;
    ASSERT_EQ(alternativeRoutes[0].route[0].start, edges[0].start);
    ASSERT_EQ(alternativeRoutes[0].route[0].end, edges[0].end);
    ASSERT_EQ(alternativeRoutes[0].route[1].start, edges[1].start);
    ASSERT_EQ(alternativeRoutes[0].route[1].end, edges[1].end);
    ASSERT_EQ(alternativeRoutes[0].route[2].start, edges[2].start);
    ASSERT_EQ(alternativeRoutes[0].route[2].end, edges[2].end);
    ASSERT_EQ(alternativeRoutes.size(), 1);
}

// Tests for initializeCars method
TEST_F(TrafficGraphTest, InitializeCars) {
    unsigned int numCars = 5;
    double minDistanceThreshold = 1;
    std::vector<Car> cars;

    graph.initializeCars(cars, numCars, minDistanceThreshold);
    ASSERT_EQ(cars.size(), numCars);

    for (size_t i = 0; i < cars.size(); ++i) {
        ASSERT_NE(cars[i].source, nullptr);
        ASSERT_NE(cars[i].destination, nullptr);
        ASSERT_TRUE(graph.isDistanceSufficient(*cars[i].source, *cars[i].destination, minDistanceThreshold));
    }
}

// // Tests for setCarRoute method
TEST_F(TrafficGraphTest, SetCarRoute) {
    std::vector<Car> cars;
    unsigned int numCars = 3;
    double minDistanceThreshold = 1;
    graph.initializeCars(cars, numCars, minDistanceThreshold);

    Route testRoute;
    testRoute.pathLen = 2;
    testRoute.route[0] = edges[0];
    testRoute.route[1] = edges[1];

    graph.setCarRoute(cars[0], testRoute, 0);
    ASSERT_EQ(cars[0].possibleRoutes[0].pathLen, testRoute.pathLen);
    ASSERT_EQ(&(cars[0].possibleRoutes[0].route[0]), &(testRoute.route[0]));
    ASSERT_EQ(&(cars[0].possibleRoutes[0].route[1]), &(testRoute.route[1]));
}

// // Tests for updateCarRoute method
TEST_F(TrafficGraphTest, UpdateCarRoute) {
    std::vector<Car> cars;
    unsigned int numCars = 3;
    double minDistanceThreshold = 1;
    graph.initializeCars(cars, numCars, minDistanceThreshold);

    Route newRoute;
    newRoute.pathLen = 3;
    newRoute.route[0] = edges[0];
    newRoute.route[1] = edges[1];
    newRoute.route[2] = edges[2];

    graph.setCarRoute(cars[0], newRoute, 0);
    graph.updateCarRoute(cars[0], newRoute, 0);

    ASSERT_EQ(cars[0].possibleRoutes[0].pathLen, newRoute.pathLen);

    printRoute(cars[0].possibleRoutes[0]);
    printRoute(newRoute);

    ASSERT_EQ(&(cars[0].possibleRoutes[0].route[0]), &(newRoute.route[0]));
    ASSERT_EQ(&(cars[0].possibleRoutes[0].route[1]), &(newRoute.route[1]));
    ASSERT_EQ(&(cars[0].possibleRoutes[0].route[2]), &(newRoute.route[2]));
}

// // Tests for calculateRouteCost method
TEST_F(TrafficGraphTest, CalculateRouteCost) {
    Route route;
    route.pathLen = 2;
    route.route[0] = edges[0];
    route.route[1] = edges[1];

    double expectedCost = edges[0].distance + edges[1].distance;
    double actualCost = graph.calculateRouteCost(route);
    ASSERT_NEAR(actualCost, expectedCost, 1e-6); // Use a small tolerance for floating-point comparison
}

// // Tests for calculateSharedCost method
TEST_F(TrafficGraphTest, CalculateSharedCost) {
    Route route1, route2;
    route1.pathLen = route2.pathLen = 2;
    route1.route[0] = route2.route[0] = edges[0];
    route1.route[1] = edges[1];
    route2.route[1] = edges[1];

    double expectedSharedCost = 10.0; // Assuming fixed shared cost
    double actualSharedCost = graph.calculateSharedCost(route1, route2);
    ASSERT_EQ(actualSharedCost, expectedSharedCost);
}

// // Tests for routesShareSegment method
// // Tests for routesShareSegment method
TEST_F(TrafficGraphTest, RoutesShareSegment) {
    Route route1, route2;
    route1.pathLen = route2.pathLen = 2;
    route1.route[0] = route2.route[0] = edges[0];
    route1.route[1] = edges[1];
    route2.route[1] = edges[2];

    ASSERT_TRUE(graph.routesShareSegment(route1, route2));
    route2.route[0] = edges[2]; // Change one edge to ensure routes don't share a segment
    ASSERT_FALSE(graph.routesShareSegment(route1, route2));
}
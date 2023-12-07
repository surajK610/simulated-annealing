#ifndef TRAFFIC_GRAPH_H
#define TRAFFIC_GRAPH_H

#include <vector>
#include <set>
#include <map>
#include <utility> // for std::pair
#include <unordered_set>

struct Point {
    short unsigned int x, y;

    bool operator<(const Point& other) const {
        return (x < other.x) || (x == other.x && y < other.y);
    }

    bool operator==(const Point& other) const {
        return (x == other.x) && (y == other.y);
    }
};

struct Edge {
    Point* start;
    Point* end;
    double distance;

    bool operator<(const Edge& other) const {
        return distance < other.distance;
    }
     bool operator==(const Edge& other) const {
        return ((start == other.start && end == other.end) || (start == other.end && end == other.start) && (distance == other.distance));
    }
};

struct EdgeHash {
    size_t operator()(const Edge& edge) const {
        size_t h1 = std::hash<Point*>()(edge.start);
        size_t h2 = std::hash<Point*>()(edge.end);
        size_t h3 = std::hash<double>()(edge.distance);

        return h1 ^ h2 ^ h3;
    }
};

struct Route {
    static const short unsigned int MAX_PATH_LEN = 300;
    short unsigned int pathLen;
    Edge route[MAX_PATH_LEN]; // ideally want pointers stored in this
};

struct Car {
    static const short unsigned int MAX_POSSIBLE_ROUTES = 300;
    short unsigned int id;
    Point* source;
    Point* destination;
    Route possibleRoutes[MAX_POSSIBLE_ROUTES]; // ideally want pointers stored in this
};

class TrafficGraph {
  
private:
    std::vector<Point>* points;
    std::vector<Edge>* edges;

    double calculateDistance(const Point& a, const Point& b);
    void addClosestEdges(std::set<Edge>& mstEdges);
    void addPoint(const Point& point);
    void addEdge(Point* start, Point* end);
    void findAllPathsUtil(Point* current, Point* destination, std::vector<Edge>& path, std::vector<Route>& allPaths, std::unordered_set<Point*>& visited);

public:
    TrafficGraph() : points(nullptr), edges(nullptr) {}
    TrafficGraph(std::vector<Point>* initPoints, std::vector<Edge>* initEdges)
        : points(initPoints), edges(initEdges) {}
    void add_random_edges(const std::set<Edge>& mstEdges, unsigned int additionalEdges);
    void initializePoints(unsigned int numPoints, unsigned int xBound, unsigned int yBound); // new method declaration
    void initializeGraph(unsigned int numPoints, unsigned int additionalEdges, unsigned int xBound, unsigned int yBound);
    std::vector<Route> findAllPaths(Point* source, Point* destination);
    std::vector<Route> findAlternativePaths(Point* source, Point* destination);
    
    Route* findShortestPath(Point* source, Point* destination);
    std::vector<Point>* getPoints() { return points; }
    std::vector<Edge>* getEdges() { return edges; }

    // Car methods
    void initializeCars(std::vector<Car>& cars, unsigned int numCars, double minDistanceThreshold);
    void setCarRoute(Car& car, const Route& route, unsigned int routeIndex);
    void updateCarRoute(Car& car, const Route& newRoute, unsigned int routeIndex);

    // Helper methods
    Point* getRandomPoint();
    bool isDistanceSufficient(const Point& a, const Point& b, double threshold);
    double** routes_to_q(const std::vector<Car>& cars);
    int calculateQMatrixSize(const std::vector<Car>& cars);
    double calculateRouteCost(const Route& route);
    double calculateSharedCost(const Route& route1, const Route& route2);
    bool routesShareSegment(const Route& route1, const Route& route2);


};

#endif

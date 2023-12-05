#ifndef TRAFFIC_GRAPH_H
#define TRAFFIC_GRAPH_H

#include <vector>
#include <set>
#include <map>
#include <utility> // for std::pair

struct Point {
    short unsigned int x, y;

    bool operator<(const Point& other) const {
        return (x < other.x) || (x == other.x && y < other.y);
    }
};

struct Edge {
    Point* start;
    Point* end;
    double distance;

    bool operator<(const Edge& other) const {
        return distance < other.distance;
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
    std::vector<Point> points;
    std::vector<Edge> edges;

    // Helpers for Prim's algorithm
    double calculateDistance(const Point& a, const Point& b);
    void addClosestEdges(std::set<Edge>& mstEdges);
    void initializePoints(unsigned int numPoints, unsigned int xBound, unsigned int yBound); // new method declaration

public:
    TrafficGraph() = default;
    void initializeGraph(unsigned int numPoints, unsigned int additionalEdges, unsigned int xBound, unsigned int yBound); // modified method signature
    void addPoint(const Point& point);
    void addEdge(Point* start, Point* end);
    std::vector<Route> findAllPaths(Point* source, Point* destination);
    std::vector<Route> findAlternativePaths(const std::vector<Route>& allPaths);
    Route* findShortestPath(Point* source, Point* destination);
};

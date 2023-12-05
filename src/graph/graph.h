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

    double calculateDistance(const Point& a, const Point& b);
    void addClosestEdges(std::set<Edge>& mstEdges);
    void initializePoints(unsigned int numPoints, unsigned int xBound, unsigned int yBound);
    void findAllPathsUtil(Point* current, Point* destination, std::vector<Edge>& path, std::vector<Route>& allPaths, std::unordered_set<Point*>& visited);

public:
    TrafficGraph() = default;

    void TrafficGraph::initializeGraph(unsigned int numPoints, unsigned int additionalEdges, unsigned int xBound, unsigned int yBound);
    void TrafficGraph::addPoint(const Point& point);
    void TrafficGraph::addEdge(Point* start, Point* end);

    std::vector<Route> TrafficGraph::findAllPaths(Point* source, Point* destination);
    std::vector<Route> TrafficGraph::findAlternativePaths(Point* source, Point* destination);
    Route* TrafficGraph::findShortestPath(Point* source, Point* destination);
};

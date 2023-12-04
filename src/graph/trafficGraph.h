#ifndef TRAFFIC_GRAPH_H
#define TRAFFIC_GRAPH_H

struct Point {
    short unsigned int x, y;
};

struct Edge {
    Point* start;
		Point* end;
		double distance;
// Prim's Algorithm (MST) (x <= 25c2)
// Randomly Sample some edges (after that) 
};

struct Route {
	static const short unsigned int MAX_ROUTES = 300;
	short unsigned int pathLen;
	Edge route[MAX_ROUTES]; // ideally want pointers stored in this
};

struct Car {
    int id;
    Point* source;
    Point* destination;
    std::vector<Route> possibleRoutes;
};

class TrafficGraph {
    std::vector<Point> nodes;
    std::vector<Edge> edges;

    void mapEndpoints(Car& car);
    std::vector<Route> findAllPaths(Point* source, Point* destination);
    std::vector<Route> findAlternativePaths(const std::vector<Route>& allPaths);
};
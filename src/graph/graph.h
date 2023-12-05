//trafficGraph.h

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
	static const short unsigned int MAX_PATH_LEN = 300;
	short unsigned int pathLen;
	Edge route[MAX_PATH_LEN]; // ideally want pointers stored in this
};

struct Car {
    static const short unsigned int MAX_POSSIBLE_ROUTES= 300;
    short unsigned int id;
    Point* source;
    Point* destination;
    Route possibleRoutes[MAX_POSSIBLE_ROUTES]; // ideally want pointers stored in this
};

class TrafficGraph {
    
    std::vector<Point> points;
    std::vector<Edge> edges;

    std::vector<Route> findAllPaths(Point* source, Point* destination);
    std::vector<Route> findAlternativePaths(const std::vector<Route>& allPaths);
    Route* findShortestPath(Point* source, Point* destination);
};
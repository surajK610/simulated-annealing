#include "trafficGraph.h"
#include <random>
#include <limits>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <random>

double TrafficGraph::calculateDistance(const Point& a, const Point& b) {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}

void TrafficGraph::addClosestEdges(std::set<Edge>& mstEdges) {
    std::unordered_set<Point*> inMst;
    for (const auto& edge : mstEdges) {
        inMst.insert(edge.start);
        inMst.insert(edge.end);
    }

    std::vector<Edge> closestEdges;
    for (const auto& point : points) {
        if (inMst.find(&point) == inMst.end()) continue;

        double minDistance = std::numeric_limits<double>::max();
        Point* closestPoint = nullptr;

        for (const auto& other : points) {
            if (&point == &other || inMst.find(&other) != inMst.end()) continue;
            double distance = calculateDistance(point, other);
            if (distance < minDistance) {
                minDistance = distance;
                closestPoint = const_cast<Point*>(&other);
            }
        }

        if (closestPoint) {
            closestEdges.emplace_back(Edge{const_cast<Point*>(&point), closestPoint, minDistance});
        }
    }

    // Randomly add a subset of these edges
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(closestEdges.begin(), closestEdges.end(), g);

    // Assuming we want to add a fixed number or a percentage of these edges
    size_t numEdgesToAdd = closestEdges.size() / 10; // for example, 10%
    for (size_t i = 0; i < numEdgesToAdd; ++i) {
        edges.push_back(closestEdges[i]);
    }
}

void TrafficGraph::initializePoints(unsigned int numPoints, unsigned int xBound, unsigned int yBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> disX(0, xBound);
    std::uniform_int_distribution<> disY(0, yBound);
    std::unordered_set<std::string> uniqueCheck;

    while (points.size() < numPoints) {
        Point newPoint;
        newPoint.x = disX(gen);
        newPoint.y = disY(gen);

        std::string pointKey = std::to_string(newPoint.x) + "-" + std::to_string(newPoint.y);
        if (uniqueCheck.find(pointKey) == uniqueCheck.end()) {
            points.push_back(newPoint);
            uniqueCheck.insert(pointKey);
            std::cout << "Generated unique point: (" << newPoint.x << ", " << newPoint.y << ")\n";
        }
    }
}

void TrafficGraph::initializeGraph(unsigned int numPoints, unsigned int additionalEdges, unsigned int xBound, unsigned int yBound) {
    std::cout << "Initializing points...\n";
    initializePoints(numPoints, xBound, yBound);

    if (points.size() < 2) {
        std::cout << "Insufficient points to form a graph.\n";
        return;
    }

    std::cout << "Initializing graph using Prim's algorithm...\n";

    // Minimum Spanning Tree using Prim's Algorithm
    auto comp = [](const Edge& e1, const Edge& e2) { return e1.distance > e2.distance; };
    std::priority_queue<Edge, std::vector<Edge>, decltype(comp)> pq(comp);

    // Start from the first point
    std::unordered_set<Point*> inMST;
    inMST.insert(&points[0]);

    // Add all edges from the starting point to the priority queue
    for (size_t i = 1; i < points.size(); ++i) {
        double dist = calculateDistance(points[0], points[i]);
        pq.push(Edge{&points[0], &points[i], dist});
    }

    std::set<Edge> mstEdges;

    while (!pq.empty() && inMST.size() < points.size()) {
        Edge smallestEdge = pq.top();
        pq.pop();

        // Check if the edge forms a cycle
        if (inMST.find(smallestEdge.end) != inMST.end()) continue;

        // Add edge to MST
        mstEdges.insert(smallestEdge);
        inMST.insert(smallestEdge.end);
        std::cout << "Added edge between (" << smallestEdge.start->x << "," << smallestEdge.start->y << ") and (" << smallestEdge.end->x << "," << smallestEdge.end->y << ") with distance " << smallestEdge.distance << ".\n";

        // Add new edges to the priority queue
        for (const auto& point : points) {
            if (inMST.find(&point) == inMST.end()) {
                double dist = calculateDistance(*smallestEdge.end, point);
                pq.push(Edge{smallestEdge.end, const_cast<Point*>(&point), dist});
            }
        }
    }

    std::cout << "Prim's algorithm completed. MST formed.\n";

    // Add additional edges
    addClosestEdges(mstEdges);
    std::cout << "Additional edges added.\n";
}

void TrafficGraph::addPoint(const Point& point) {
    points.push_back(point);
}

void TrafficGraph::addEdge(Point* start, Point* end) {
    double distance = calculateDistance(*start, *end);
    edges.push_back(Edge{start, end, distance});
}


std::vector<Route> TrafficGraph::findAllPaths(Point* source, Point* destination);
std::vector<Route> TrafficGraph::findAlternativePaths(Point* source, Point* destination);


Route* TrafficGraph::findShortestPath(Point* source, Point* destination) {
    std::unordered_map<Point*, double> distances;
    
    std::unordered_map<Point*, Edge*> previous;

    for (auto& point : points) {
        distances[&point] = std::numeric_limits<double>::infinity();
    }

    auto cmp = [&distances](Point* left, Point* right) { return distances[left] > distances[right]; };
    std::priority_queue<Point*, std::vector<Point*>, decltype(cmp)> queue(cmp);

    distances[source] = 0;
    queue.push(source);

    while (!queue.empty()) {
        Point* current = queue.top();
        queue.pop();

        // If we reached the destination, reconstruct and return the route
        if (current == destination) {
            return reconstructRoute(source, destination, previous);
        }

        // Relaxation step for each edge
        for (auto& edge : edges) {
            if (edge.start == current) {
                Point* neighbor = edge.end;
                double newDist = distances[current] + edge.distance;
                if (newDist < distances[neighbor]) {
                    distances[neighbor] = newDist;
                    previous[neighbor] = &edge;
                    queue.push(neighbor);
                }
            }
        }
    }

    // Destination is not reachable
    return nullptr;
}

Route* reconstructRoute(Point* source, Point* destination, std::unordered_map<Point*, Edge*>& previous) {
    // Construct the shortest route from destination to source using previous map
    std::vector<Edge*> path;
    for (Point* at = destination; at != source; at = previous[at]->start) {
        path.push_back(previous[at]);
    }
    std::reverse(path.begin(), path.end());

    // Convert path to Route and return it
    Route* shortestRoute = new Route();
    shortestRoute->pathLen = path.size();
    for (size_t i = 0; i < path.size(); ++i) {
        shortestRoute->route[i] = *path[i];
    }

    return shortestRoute;
}

#include "graph.h"
#include <random>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>

#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>

double TrafficGraph::calculateDistance(const Point& a, const Point& b) {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}

void TrafficGraph::addClosestEdges(std::set<Edge>& mstEdges) {
    std::unordered_set<Point*> inMST;
    for (const auto& edge : mstEdges) {
        inMST.insert(edge.start);
        inMST.insert(edge.end);
    }

    std::vector<Edge> closestEdges;
    for (auto& point : points) {
        Point* pointPtr = &point;
        if (inMST.find(pointPtr) != inMST.end()) continue;

        double minDistance = std::numeric_limits<double>::max();
        Point* closestPoint = nullptr;

        for (auto& other : points) {
            Point* otherPtr = &other;
            if (pointPtr == otherPtr || inMST.find(otherPtr) != inMST.end()) continue;
            double distance = calculateDistance(*pointPtr, *otherPtr);
            if (distance < minDistance) {
                minDistance = distance;
                closestPoint = otherPtr;
            }
        }

        if (closestPoint) {
            closestEdges.emplace_back(Edge{pointPtr, closestPoint, minDistance});
        }
    }

    // Randomly add a subset of these edges
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(closestEdges.begin(), closestEdges.end(), g);

    size_t numEdgesToAdd = closestEdges.size() / 10; // Example: 10%
    for (size_t i = 0; i < numEdgesToAdd; ++i) {
        edges.push_back(closestEdges[i]);
    }

    std::cerr << "Adding the Additional Edges To The Edge List (ClostestEdgeFunction)...\n";

    // Log the edges that were added
    for (size_t i = 0; i < numEdgesToAdd; ++i) {
        const Edge& addedEdge = closestEdges[i];
        std::cerr << i + 1 << ": Added edge between (" << addedEdge.start->x << "," << addedEdge.start->y
                  << ") and (" << addedEdge.end->x << "," << addedEdge.end->y << ") with distance "
                  << addedEdge.distance << ".\n";
    }
}


void TrafficGraph::initializePoints(unsigned int numPoints, unsigned int xBound, unsigned int yBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> disX(0, xBound);
    std::uniform_int_distribution<> disY(0, yBound);
    std::unordered_set<std::string> uniqueCheck;
    std::cerr << "Starting Unique Point Generation...\n";

    while (points.size() < numPoints) {
        short unsigned int x = static_cast<short unsigned int>(disX(gen));
        short unsigned int y = static_cast<short unsigned int>(disY(gen));
        Point newPoint{x, y};

        std::string pointKey = std::to_string(newPoint.x) + "-" + std::to_string(newPoint.y);
        if (uniqueCheck.insert(pointKey).second) {
            points.push_back(newPoint);
            std::cerr << "Generated unique point: (" << newPoint.x << ", " << newPoint.y << ")\n";
        }
    }
}



void TrafficGraph::initializeGraph(unsigned int numPoints, unsigned int additionalEdges, unsigned int xBound, unsigned int yBound) {
    
    // Clear the points and edges vectors to start with empty containers
    points.clear();
    edges.clear();
    std::cerr << "Points vector cleared. Size after: " << points.size() << std::endl;
    std::cerr << "Edges vector cleared. Size after: " << edges.size() << std::endl;

    std::cerr << "Initializing points...\n";
    std::unordered_set<Point*> inMST;

    initializePoints(numPoints, xBound, yBound);

    if (points.size() < 2) {
        std::cerr << "Insufficient points to form a graph.\n";
        return;
    }

    std::cerr << "Initializing graph using Prim's algorithm...\n";

    auto comp = [](const Edge& e1, const Edge& e2) { return e1.distance > e2.distance; };
    std::priority_queue<Edge, std::vector<Edge>, decltype(comp)> pq(comp);

    inMST.insert(&points[0]);

    for (size_t i = 1; i < points.size(); ++i) {
        double dist = calculateDistance(points[0], points[i]);
        pq.push(Edge{&points[0], &points[i], dist});
    }

    std::set<Edge> mstEdges;

    while (!pq.empty() && inMST.size() < points.size()) {
        Edge smallestEdge = pq.top();
        pq.pop();

        if (inMST.find(smallestEdge.end) != inMST.end()) continue;

        mstEdges.insert(smallestEdge);
        inMST.insert(smallestEdge.end);
        std::cerr << "Added edge between (" << smallestEdge.start->x << "," << smallestEdge.start->y << ") and (" << smallestEdge.end->x << "," << smallestEdge.end->y << ") with distance " << smallestEdge.distance << ".\n";

        for (auto& point : points) {
            if (inMST.find(&point) == inMST.end()) {
                double dist = calculateDistance(*smallestEdge.end, point);
                pq.push(Edge{smallestEdge.end, const_cast<Point*>(&point), dist});
            }
        }
    }

    std::cout << "Prim's algorithm completed. MST formed.\n";

    addClosestEdges(mstEdges);
    std::cout << "Additional edges added.\n";

    // Add each edge from mstEdges to the edges list
    for (const auto& edge : mstEdges) {
        addEdge(edge.start, edge.end);
    }

    std::cout << "Additional edges added.\n";
}

void TrafficGraph::addPoint(const Point& point) {
    points.push_back(point);
}

void TrafficGraph::addEdge(Point* start, Point* end) {
    double distance = calculateDistance(*start, *end);
    edges.push_back(Edge{start, end, distance});
}

void  TrafficGraph::findAllPathsUtil(Point* current, Point* destination, std::vector<Edge>& path, std::vector<Route>& allPaths, std::unordered_set<Point*>& visited) {
    if (current == destination) {
        Route route;
        route.pathLen = path.size();
        for (size_t i = 0; i < path.size(); ++i) {
            route.route[i] = path[i];
        }
        allPaths.push_back(route);
        return;
    }

    visited.insert(current);

    for (auto& edge : edges) {
        if (edge.start == current && visited.find(edge.end) == visited.end()) {
            path.push_back(edge);
            findAllPathsUtil(edge.end, destination, path, allPaths, visited);
            path.pop_back();
        } else if (edge.end == current && visited.find(edge.start) == visited.end()) {
            path.push_back(edge);
            findAllPathsUtil(edge.start, destination, path, allPaths, visited);
            path.pop_back();
        }
    }

    visited.erase(current);
}

std::vector<Route>  TrafficGraph::findAllPaths(Point* source, Point* destination) {
        std::vector<Route> allPaths;
        std::vector<Edge> path;
        std::unordered_set<Point*> visited;

        findAllPathsUtil(source, destination, path, allPaths, visited);
        return allPaths;
}

double jaccardSimilarity(Route& route1, Route& route2) {
        std::unordered_set<Edge*> edgesInRoute1;
        std::unordered_set<Edge*> edgesInRoute2;

        for (int i = 0; i < route1.pathLen; ++i) {
            edgesInRoute1.insert(&route1.route[i]);
        }

        for (int i = 0; i < route2.pathLen; ++i) {
            edgesInRoute2.insert(&route2.route[i]);
        }

        int intersectionSize = 0;
        for (auto edge : edgesInRoute1) {
            if (edgesInRoute2.find(edge) != edgesInRoute2.end()) {
                intersectionSize++;
            }
        }

        int unionSize = edgesInRoute1.size() + edgesInRoute2.size() - intersectionSize;
        return static_cast<double>(intersectionSize) / unionSize;
}

Route* reconstructRoute(Point* source, Point* destination, std::unordered_map<Point*, Edge*>& previous) {
    std::vector<Edge*> path;
    for (Point* at = destination; at != source; at = previous[at]->start) {
        path.push_back(previous[at]);
    }
    std::reverse(path.begin(), path.end());

    Route* shortestRoute = new Route();
    shortestRoute->pathLen = path.size();
    for (size_t i = 0; i < path.size(); ++i) {
        shortestRoute->route[i] = *path[i];
    }

    return shortestRoute;
}

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

        if (current == destination) {
            std::cerr << "Found shortest path from (" << source->x << "," << source->y << ") to (" << destination->x << "," << destination->y << ") with distance " << distances[current] << ".\n";
            return reconstructRoute(source, destination, previous);
        }
        std::cerr << "Made here\n";

        for (auto& edge : edges) {
            Point* neighbor = nullptr;
            if (edge.start == current) {
                neighbor = edge.end;
            } else if (edge.end == current) {
                neighbor = edge.start; // Consider the edge in the reverse direction for undirected graph
            }

            if (neighbor) {
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

std::vector<Route> TrafficGraph::findAlternativePaths(Point* source, Point* destination) {
    Route* shortestRoute = findShortestPath(source, destination);
    std::vector<Route> allRoutes = findAllPaths(source, destination);

    std::vector<std::pair<double, Route>> similarityScores;

    for (auto& route : allRoutes) {
        double similarity = jaccardSimilarity(*shortestRoute, route);
        similarityScores.push_back(std::make_pair(similarity, route));
    }

    std::sort(similarityScores.begin(), similarityScores.end(), [](const std::pair<double, Route>& a, const std::pair<double, Route>& b) {
        return a.first < b.first;
    });

    std::vector<Route> alternativeRoutes;
    for (size_t i = 0; i < similarityScores.size() && alternativeRoutes.size() < 2; ++i) {
        if (similarityScores[i].first != 0) {
            alternativeRoutes.push_back(similarityScores[i].second);
        }
    }

    delete shortestRoute;
    return alternativeRoutes;
}
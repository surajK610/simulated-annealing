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

    auto comp = [](const Edge& e1, const Edge& e2) { return e1.distance > e2.distance; };
    std::priority_queue<Edge, std::vector<Edge>, decltype(comp)> pq(comp);

    std::unordered_set<Point*> inMST;
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
        std::cout << "Added edge between (" << smallestEdge.start->x << "," << smallestEdge.start->y << ") and (" << smallestEdge.end->x << "," << smallestEdge.end->y << ") with distance " << smallestEdge.distance << ".\n";

        for (const auto& point : points) {
            if (inMST.find(&point) == inMST.end()) {
                double dist = calculateDistance(*smallestEdge.end, point);
                pq.push(Edge{smallestEdge.end, const_cast<Point*>(&point), dist});
            }
        }
    }

    std::cout << "Prim's algorithm completed. MST formed.\n";

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

void TrafficGraph::findAllPathsUtil(Point* current, Point* destination, std::vector<Edge>& path, std::vector<Route>& allPaths, std::unordered_set<Point*>& visited) {
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
        // Check if the edge is connected to the current point and the other point is not visited
        if (edge.start == current && visited.find(edge.end) == visited.end()) {
            path.push_back(edge);
            findAllPathsUtil(edge.end, destination, path, allPaths, visited);
            path.pop_back();
        } else if (edge.end == current && visited.find(edge.start) == visited.end()) {
            // Since the edges are undirected, also check the other direction
            path.push_back(edge);
            findAllPathsUtil(edge.start, destination, path, allPaths, visited);
            path.pop_back();
        }
    }

    visited.erase(current);
}

std::vector<Route> findAllPaths(Point* source, Point* destination) {
        std::vector<Route> allPaths;
        std::vector<Edge> path;
        std::unordered_set<Point*> visited;

        findAllPathsUtil(source, destination, path, allPaths, visited);
        return allPaths;
}

double jaccardSimilarity(const Route& route1, const Route& route2) {
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

std::vector<Route> findAlternativePaths(Point* source, Point* destination) {
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
            return reconstructRoute(source, destination, previous);
        }

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
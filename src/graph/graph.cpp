//trafficGraph.cpp

#include "graph.h"
#include <random>
#include <limits>
#include <omp.h>

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
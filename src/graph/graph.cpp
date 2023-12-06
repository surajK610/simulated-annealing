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
    for (auto& point : *points) {
        Point* pointPtr = &point;
        if (inMST.find(pointPtr) != inMST.end()) continue;

        double minDistance = std::numeric_limits<double>::max();
        Point* closestPoint = nullptr;

        for (auto& other : *points) {
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
        edges->push_back(closestEdges[i]);
    }

    std::cerr << "Adding the Additional Edges To The Edge List (ClostestEdgeFunction)...\n";
    std::cerr << "Adding" << numEdgesToAdd << " edges to the edge list.\n";

    // Log the edges that were added
    for (size_t i = 0; i < numEdgesToAdd; ++i) {
        const Edge& addedEdge = closestEdges[i];
        std::cerr << i + 1 << ": Added edge between (" << addedEdge.start->x << "," << addedEdge.start->y
                  << ") and (" << addedEdge.end->x << "," << addedEdge.end->y << ") with distance "
                  << addedEdge.distance << ".\n";
    }
}


void TrafficGraph::add_random_edges(const std::set<Edge>& mstEdges, unsigned int additionalEdges) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, points->size() - 1);

    unsigned int addedEdges = 0, numTries = 0;
    const double MIN_DISTANCE_THRESHOLD = 30.0;  // Adjust this threshold as needed
    const unsigned int MAX_TRIES = 200;

    std::cerr << "Starting to add random edges...\n";
    while (addedEdges < additionalEdges) {
        if (numTries >= MAX_TRIES) {
            std::cerr << "Exceeded maximum number of tries, stopping...\n";
            break;
        }

        int idx1 = dis(gen);
        int idx2 = dis(gen);

        //std::cerr << "Random indices selected: idx1 = " << idx1 << ", idx2 = " << idx2 << "\n";

        // Ensure we don't select the same point
        if (idx1 == idx2) {
            numTries++;
            //std::cerr << "Selected the same point for both ends of edge, retrying...\n";
            continue;
        }

        Point* start = &((*points)[idx1]);
        Point* end = &((*points)[idx2]);
        double distance = calculateDistance(*start, *end);

        //std::cerr << "Calculated distance: " << distance << "\n";

        if (distance > MIN_DISTANCE_THRESHOLD) {
            //std::cerr << "Distance below threshold, skipping this pair...\n";
            numTries++;
            continue;
        }

        // Check if this edge already exists in mstEdges
        Edge newEdge = {start, end, distance};
        auto it = std::find_if(mstEdges.begin(), mstEdges.end(), [&newEdge](const Edge& e) {
            return (*e.start == *newEdge.start && *e.end == *newEdge.end) || (*e.start == *newEdge.end && *e.end == *newEdge.start);
        });

        if (it == mstEdges.end()) {
            // Edge not in mstEdges, so add it
            addEdge(start, end);
            std::cerr << "Random edge added between (" << start->x << "," << start->y << ") and (" << end->x << "," << end->y << ") with distance " << distance << ".\n";
            ++addedEdges;
        } else {
            numTries++;
            std::cerr << "Edge already exists in MST, skipping...\n";
        }
    }

    std::cerr << "Finished adding random edges. Total added: " << addedEdges << "\n";
}


void TrafficGraph::initializePoints(unsigned int numPoints, unsigned int xBound, unsigned int yBound) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> disX(0, xBound);
    std::uniform_int_distribution<> disY(0, yBound);
    std::unordered_set<std::string> uniqueCheck;
    std::cerr << "Starting Unique Point Generation...\n";

    while (points->size() < numPoints) {
        short unsigned int x = static_cast<short unsigned int>(disX(gen));
        short unsigned int y = static_cast<short unsigned int>(disY(gen));
        Point newPoint{x, y};

        std::string pointKey = std::to_string(newPoint.x) + "-" + std::to_string(newPoint.y);
        if (uniqueCheck.insert(pointKey).second) {
            points->push_back(newPoint);
            // std::cerr << "Generated unique point: (" << newPoint.x << ", " << newPoint.y << ")\n";
        }
    }
}



void TrafficGraph::initializeGraph(unsigned int numPoints, unsigned int additionalEdges, unsigned int xBound, unsigned int yBound) {
    
    // Clear the points and edges vectors to start with empty containers
    points->clear();
    edges->clear();
    std::cerr << "Points vector cleared. Size after: " << points->size() << std::endl;
    std::cerr << "Edges vector cleared. Size after: " << edges->size() << std::endl;

    std::cerr << "Initializing points...\n";
    std::unordered_set<Point*> inMST;

    initializePoints(numPoints, xBound, yBound);

    if (points->size() < 2) {
        std::cerr << "Insufficient points to form a graph.\n";
        return;
    }

    std::cerr << "Initializing graph using Prim's algorithm...\n";

    auto comp = [](const Edge& e1, const Edge& e2) { return e1.distance > e2.distance; };
    std::priority_queue<Edge, std::vector<Edge>, decltype(comp)> pq(comp);

    inMST.insert(&(*points)[0]);

    for (size_t i = 1; i < points->size(); ++i) {
        double dist = calculateDistance((*points)[0], (*points)[i]);
        pq.push(Edge{&(*points)[0], &(*points)[i], dist});
    }

    std::set<Edge> mstEdges;

    while (!pq.empty() && inMST.size() < points->size()) {
        Edge smallestEdge = pq.top();
        pq.pop();

        if (inMST.find(smallestEdge.end) != inMST.end()) continue;

        mstEdges.insert(smallestEdge);
        inMST.insert(smallestEdge.end);
        // std::cerr << "Added edge between (" << smallestEdge.start->x << "," << smallestEdge.start->y << ") and (" << smallestEdge.end->x << "," << smallestEdge.end->y << ") with distance " << smallestEdge.distance << ".\n";

        for (auto& point : *points) {
            if (inMST.find(&point) == inMST.end()) {
                double dist = calculateDistance(*smallestEdge.end, point);
                pq.push(Edge{smallestEdge.end, const_cast<Point*>(&point), dist});
            }
        }
    }

    std::cerr << "Prim's algorithm completed. MST formed.\n";

    add_random_edges(mstEdges, additionalEdges);
    std::cerr << "Additional edges added.\n";

    // Add each edge from mstEdges to the edges list
    for (const auto& edge : mstEdges) {
        addEdge(edge.start, edge.end);
    }

    std::cerr << "Additional edges added.\n";
}

void TrafficGraph::addPoint(const Point& point) {
    points->push_back(point);
}

void TrafficGraph::addEdge(Point* start, Point* end) {
    double distance = calculateDistance(*start, *end);
    edges->push_back(Edge{start, end, distance});
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

    for (auto& edge : *edges) {
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
    // std::cerr << "Source" << source << " " << source->x << " " << source->y << std::endl;
    // std::cerr << "Destination" << destination << " " << destination->x << " " << destination->y << std::endl;

    for (auto& point : *points) {
        distances[&point] = std::numeric_limits<double>::infinity();
        // std::cerr << &point << " " << point.x << " " << point.y << std::endl;
    }
    
    auto cmp = [&distances](Point* left, Point* right) { return distances[left] > distances[right]; };
    std::priority_queue<Point*, std::vector<Point*>, decltype(cmp)> queue(cmp);

    // std::cerr << "Source" << source << " " << source->x << " " << source->y << std::endl;
    distances[source] = 0;
    queue.push(source);
    // std::cerr << "Made here1\n";
    while (!queue.empty()) {
        Point* current = queue.top();
        queue.pop();
        // std::cerr << current->x << " " << current->y << std::endl;

        if (current == destination) {
            // std::cerr << "Found shortest path from (" << source->x << "," << source->y << ") to (" << destination->x << "," << destination->y << ") with distance " << distances[current] << ".\n";
            return reconstructRoute(source, destination, previous);
        }
        // std::cerr << "Made here\n";

        for (auto& edge : *edges) {
            Point* neighbor = nullptr;
            if (edge.start == current) {
                neighbor = edge.end;
            } else if (edge.end == current) {
                neighbor = edge.start; // Consider the edge in the reverse direction for undirected graph
            }

            if (neighbor) {
                // std::cerr << "Neighbor: " << neighbor << " " << neighbor->x << " " << neighbor->y << "\n";
                double newDist = distances[current] + edge.distance;
                // std::cerr << "New distance: " << newDist << "\n";

                if (newDist < distances[neighbor]) {
                    distances[neighbor] = newDist;
                    previous[neighbor] = &edge;
                    queue.push(neighbor);
                    // std::cerr << "Queue size: " << queue.size() << "\n";
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

//////////////////////////////////////////////////////////////
// TrafficRoute Methods (Cars)
//////////////////////////////////////////////////////////////



double** TrafficGraph::routes_to_q(const std::vector<Car>& cars) {
    // Determine the size of the Q matrix
    int size = calculateQMatrixSize(cars); // You need to implement this
    double** Q = new double*[size];

    for (int i = 0; i < size; ++i) {
        Q[i] = new double[size](); // Initialize to zero
    }

    try {
        for (const auto& car : cars) {
            for (int i = 0; i < Car::MAX_POSSIBLE_ROUTES; ++i) {
                // Set diagonal based on route cost
                double routeCost = calculateRouteCost(car.possibleRoutes[i]); // Implement this
                Q[i][i] = routeCost;

                // Set off-diagonal based on shared segments
                for (int j = i + 1; j < Car::MAX_POSSIBLE_ROUTES; ++j) {
                    if (routesShareSegment(car.possibleRoutes[i], car.possibleRoutes[j])) { // Implement this
                        double sharedCost = calculateSharedCost(car.possibleRoutes[i], car.possibleRoutes[j]); // Implement this
                        Q[i][j] = Q[j][i] = sharedCost;
                    }
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error in routes_to_q: " << e.what() << std::endl;
        // Cleanup if needed
        for (int i = 0; i < size; ++i) {
            delete[] Q[i];
        }
        delete[] Q;
        throw; // rethrow the exception
    }

    return Q;
}


int TrafficGraph::calculateQMatrixSize(const std::vector<Car>& cars) {
    if (cars.empty()) return 0;
    int numRoutesPerCar = Car::MAX_POSSIBLE_ROUTES; // Assuming a constant number of possible routes per car
    return cars.size() * numRoutesPerCar;
}

double TrafficGraph::calculateRouteCost(const Route& route) {
    double cost = 0.0;
    for (int i = 0; i < route.pathLen - 1; ++i) {
        Point* start = route.route[i].start; // Corrected
        Point* end = route.route[i].end;     // Corrected
        cost += calculateDistance(*start, *end);
    }
    return cost;
}

double TrafficGraph::calculateSharedCost(const Route& route1, const Route& route2) {
    double sharedCost = 0.0;
    if (routesShareSegment(route1, route2)) {
        // Cost can be a constant or based on the length of shared segments
        sharedCost = 10.0; // Example fixed cost
    }
    return sharedCost;
}

bool TrafficGraph::routesShareSegment(const Route& route1, const Route& route2) {
    std::unordered_set<const Edge*> edgesInRoute1; // Type changed to const Edge*
    for (int i = 0; i < route1.pathLen; ++i) {
        edgesInRoute1.insert(&route1.route[i]);
    }

    for (int i = 0; i < route2.pathLen; ++i) {
        if (edgesInRoute1.find(&route2.route[i]) != edgesInRoute1.end()) {
            return true;
        }
    }
    return false;
}

// Implementation of getRandomPoint
Point* TrafficGraph::getRandomPoint() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, points->size() - 1);
    int index = dis(gen);
    return &((*points)[index]);
}

// Implementation of isDistanceSufficient
bool TrafficGraph::isDistanceSufficient(const Point& a, const Point& b, double threshold) {
    double distance = calculateDistance(a, b);
    return distance >= threshold;
}

// Implementation of initializeCars
void TrafficGraph::initializeCars(std::vector<Car>& cars, unsigned int numCars, double minDistanceThreshold) {
    cars.clear();
    cars.reserve(numCars);
    
    std::random_device rd;
    std::mt19937 gen(rd());

    for (unsigned int i = 0; i < numCars; ++i) {
        Car car;
        car.id = i;

        Point* src = nullptr;
        Point* dest = nullptr;
        do {
            src = getRandomPoint();
            dest = getRandomPoint();
        } while (!isDistanceSufficient(*src, *dest, minDistanceThreshold));

        car.source = src;
        car.destination = dest;
        cars.push_back(car);

        std::cerr << "Car ID " << car.id << " assigned source (" << src->x << ", " << src->y 
                  << ") and destination (" << dest->x << ", " << dest->y << ").\n";
    }
}

void TrafficGraph::setCarRoute(Car& car, const Route& route, unsigned int routeIndex) {
    if (routeIndex < Car::MAX_POSSIBLE_ROUTES) {
        car.possibleRoutes[routeIndex] = route;
    } else {
        std::cerr << "Error: Invalid route index in setCarRoute for car ID " << car.id << std::endl;
    }
}

void TrafficGraph::updateCarRoute(Car& car, const Route& newRoute, unsigned int routeIndex) {
    if (routeIndex < Car::MAX_POSSIBLE_ROUTES) {
        car.possibleRoutes[routeIndex] = newRoute;
    } else {
        std::cerr << "Error: Invalid route index in updateCarRoute for car ID " << car.id << std::endl;
    }
}

// main.cpp
#include "graph/graph.h"
#include "qubo/qubo.h"
#include <vector>
#include <iostream>

int main(int argc, char **argv) {
    // Initialize TrafficGraph
    TrafficGraph graph;
    unsigned int numPoints = 100; // Example: number of points in the graph
    unsigned int additionalEdges = 50; // Example: number of additional edges to add
    unsigned int xBound = 1000; // Example: max X-coordinate
    unsigned int yBound = 1000; // Example: max Y-coordinate
    graph.initializeGraph(numPoints, additionalEdges, xBound, yBound);

    std::cerr << "Initialized TrafficGraph" << std::endl;

    // Initialize cars and assign routes
    std::vector<Car> cars;
    unsigned int numCars = 418; // Example: number of cars
    double minDistanceThreshold = 100.0; // Example: minimum distance threshold for source and destination
    graph.initializeCars(cars, numCars, minDistanceThreshold);

    std::cerr << "Initialized cars and assigned routes" << std::endl;

    // Convert routes to QUBO problem
    double** Q = graph.routes_to_q(cars);
    int quboSize = graph.calculateQMatrixSize(cars);

    std::cerr << "Converted routes to QUBO problem" << std::endl;

    // Solve QUBO problem
    int* bestConfiguration = new int[quboSize];
    int numSamples = 1000; // Example: number of samples for Monte Carlo simulation
    monteCarloQUBOSolver(Q, numSamples, bestConfiguration, quboSize);

    std::cerr << "Solved QUBO problem" << std::endl;

    // Print the best configuration found
    std::cerr << "Best Configuration:" << std::endl;
    for (int i = 0; i < quboSize; ++i) {
        std::cerr << "Car " << i / 3 << " Route " << i % 3 << ": " << bestConfiguration[i] << std::endl;
    }

    // Cleanup
    for (int i = 0; i < quboSize; ++i) {
        delete[] Q[i];
    }
    delete[] Q;
    delete[] bestConfiguration;

    return 0;
}

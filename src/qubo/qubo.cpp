//quboSolver.cpp

#include "qubo.h"
#include <iostream>
#include <random>
#include <limits>
#include <omp.h>

// Function to calculate the objective function value for a given configuration
double calculateObjective(double** Q, int* configuration, int size) {
    double objectiveValue = 0.0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            objectiveValue += Q[i][j] * configuration[i] * configuration[j];
        }
    }
    return objectiveValue;
}

// Monte Carlo simulation to solve the QUBO problem
void monteCarloQUBOSolver(double** Q, int numSamples, int* bestConfiguration, int size) {
    double bestObjectiveValue = std::numeric_limits<double>::infinity();
    int* localBestConfiguration = new int[size];
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);
    // try {
    //     #pragma omp parallel shared(bestObjectiveValue, bestConfiguration)
    //     {
            int* threadLocalBestConfiguration = new int[size];
            double threadLocalBestObjectiveValue = std::numeric_limits<double>::infinity();

            // #pragma omp for nowait
            for (int sample = 0; sample < numSamples; ++sample) {
                int* currentConfiguration = new int[size];

                // Reset configuration to -1, indicating no route assigned
                std::fill_n(currentConfiguration, size, -1);

                for (int car = 0; car < numCars; ++car) {
                    int routeOffset = car * numRoutesPerCar;  // Calculate the offset in the QUBO matrix for this car

                    // Check for the feasibility of switching routes
                    if (canSwitchRoute(car)) {
                        // Assign one of the available routes randomly
                        int assignedRoute = dis(gen) % numRoutesPerCar;
                        currentConfiguration[routeOffset + assignedRoute] = 1;
                    } else {
                        // If the car cannot switch routes, assign it to its current route
                        int currentRoute = getCurrentRouteIndex(car);
                        currentConfiguration[routeOffset + currentRoute] = 1;
                    }
                }

                double currentObjectiveValue = calculateObjective(Q, currentConfiguration, size);

                if (currentObjectiveValue < threadLocalBestObjectiveValue) {
                    threadLocalBestObjectiveValue = currentObjectiveValue;
                    std::swap(threadLocalBestConfiguration, currentConfiguration);
                }
                delete[] currentConfiguration;
            }

            // #pragma omp critical
            // {
                if (threadLocalBestObjectiveValue < bestObjectiveValue) {
                    bestObjectiveValue = threadLocalBestObjectiveValue;
                    std::swap(bestConfiguration, threadLocalBestConfiguration);
                }
            // }

            // delete[] threadLocalBestConfiguration;
        // }
    // } catch (const std::exception& e) {
    //     std::cerr << "Error in monteCarloQUBOSolver: " << e.what() << std::endl;
    //     delete[] localBestConfiguration;
    //     throw;
    // }

    // delete[] localBestConfiguration;
}

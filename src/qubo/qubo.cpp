//quboSolver.cpp

#include "qubo.h"
#include <iostream>
#include <random>
#include <limits>
#include <omp.h>

// Function to calculate the objective function value for a given configuration
double calculateObjective(double Q[N][N], int configuration[N]) {
    double objectiveValue = 0.0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            objectiveValue += Q[i][j] * configuration[i] * configuration[j];
        }
    }
    return objectiveValue;
}

// Monte Carlo simulation to solve the QUBO problem
void monteCarloQUBOSolver(double Q[N][N], int numSamples, int bestConfiguration[N]) {
    double bestObjectiveValue = std::numeric_limits<double>::infinity();
    int threadBestConfiguration[N];

    // #pragma omp parallel shared(bestObjectiveValue, bestConfiguration)
    // {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    int localBestConfiguration[N];
    double localBestObjectiveValue = std::numeric_limits<double>::infinity();

    // #pragma omp for nowait
    for (int sample = 0; sample < numSamples; ++sample) {
        int currentConfiguration[N];
        for (int i = 0; i < N; ++i) {
            currentConfiguration[i] = dis(gen); // Randomly assign 0 or 1
        }

        double currentObjectiveValue = calculateObjective(Q, currentConfiguration);

        // If the new configuration is better, update the local best values
        if (currentObjectiveValue < localBestObjectiveValue) {
            localBestObjectiveValue = currentObjectiveValue;
            for (int i = 0; i < N; ++i) {
                localBestConfiguration[i] = currentConfiguration[i];
            }
        }
    }

    // Update the global best configuration
    // #pragma omp critical
    // {
    if (localBestObjectiveValue < bestObjectiveValue) {
        bestObjectiveValue = localBestObjectiveValue;
        for (int i = 0; i < N; ++i) {
            bestConfiguration[i] = localBestConfiguration[i];
        }
    }
        // }
    // }
}
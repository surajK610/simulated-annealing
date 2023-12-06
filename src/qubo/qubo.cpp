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

    try {
        #pragma omp parallel shared(bestObjectiveValue, bestConfiguration)
        {
            int* threadLocalBestConfiguration = new int[size];
            double threadLocalBestObjectiveValue = std::numeric_limits<double>::infinity();

            #pragma omp for nowait
            for (int sample = 0; sample < numSamples; ++sample) {
                int* currentConfiguration = new int[size];
                for (int i = 0; i < size; ++i) {
                    currentConfiguration[i] = dis(gen); // Randomly assign 0 or 1
                }

                double currentObjectiveValue = calculateObjective(Q, currentConfiguration, size);

                if (currentObjectiveValue < threadLocalBestObjectiveValue) {
                    threadLocalBestObjectiveValue = currentObjectiveValue;
                    std::swap(threadLocalBestConfiguration, currentConfiguration);
                }

                delete[] currentConfiguration;
            }

            #pragma omp critical
            {
                if (threadLocalBestObjectiveValue < bestObjectiveValue) {
                    bestObjectiveValue = threadLocalBestObjectiveValue;
                    std::swap(bestConfiguration, threadLocalBestConfiguration);
                }
            }

            delete[] threadLocalBestConfiguration;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error in monteCarloQUBOSolver: " << e.what() << std::endl;
        delete[] localBestConfiguration;
        throw;
    }

    delete[] localBestConfiguration;
}

#include "gtest/gtest.h"
#include "../src/qubo/qubo.h" // Update the relative path to your quboSolver.h

// Helper function to create a dynamic 2D array from a static one
double** createDynamic2DArray(double staticArray[][3], int N) {
    double** dynamicArray = new double*[N];
    for (int i = 0; i < N; ++i) {
        dynamicArray[i] = new double[N];
        for (int j = 0; j < N; ++j) {
            dynamicArray[i][j] = staticArray[i][j];
        }
    }
    return dynamicArray;
}

// Helper function to delete a dynamic 2D array
void deleteDynamic2DArray(double** array, int N) {
    for (int i = 0; i < N; ++i) {
        delete[] array[i];
    }
    delete[] array;
}

// Test for CalculateObjective function
TEST(Qubo, CalculateObjective) {
    double staticQ[3][3] = {
        {1, -2, 1},
        {-2, 4, -2},
        {1, -2, 1}
    };
    double** Q = createDynamic2DArray(staticQ, 3);
    int configuration[3] = {1, 0, 1};

    double expectedValue = 4; // Calculate this by hand for your test case
    double objectiveValue = QUBO::calculateObjective(Q, configuration, 3);

    deleteDynamic2DArray(Q, 3); // Clean up dynamic array
    EXPECT_DOUBLE_EQ(expectedValue, objectiveValue);
}

// A test for monteCarloQUBOSolver function
TEST(Qubo, MonteCarloQUBOSolver) {
    double staticQ[3][3] = {
        {1, -2, 1},
        {-2, 4, -2},
        {1, -2, 1}
    };
    double** Q = createDynamic2DArray(staticQ, 3);
    int numSamples = 1000; // Number of samples for the Monte Carlo simulation
    int bestConfiguration[3];

    QUBO::monteCarloQUBOSolver(Q, 3, numSamples, bestConfiguration, 3);

    // Ensure the best configuration is a valid binary string
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(bestConfiguration[i] == 0 || bestConfiguration[i] == 1);
    }

    deleteDynamic2DArray(Q, 3); // Clean up dynamic array
}

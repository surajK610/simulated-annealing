// tests/qubo_test.cpp

#include "gtest/gtest.h"
#include "../src/qubo/qubo.h" // Update the relative path to your quboSolver.h


TEST(Qubo, CalculateObjective) {
    double Q[N][N] = {
        {1, -2, 1},
        {-2, 4, -2},
        {1, -2, 1}
    };
    int configuration[N] = {1, 0, 1};

    double expectedValue = 4; // Calculate this by hand for your test case
    double objectiveValue = calculateObjective(Q, configuration);
    EXPECT_DOUBLE_EQ(expectedValue, objectiveValue);
}

// A test for monteCarloQUBOSolver function
// This test ensures that the function returns a valid configuration
TEST(Qubo, MonteCarloQUBOSolver) {
    double Q[N][N] = {
        {1, -2, 1},
        {-2, 4, -2},
        {1, -2, 1}
    };
    int numSamples = 1000; // Number of samples for the Monte Carlo simulation
    int bestConfiguration[N];

    monteCarloQUBOSolver(Q, numSamples, bestConfiguration);

    // Ensure the best configuration is a valid binary string
    for (int i = 0; i < N; ++i) {
        EXPECT_TRUE(bestConfiguration[i] == 0 || bestConfiguration[i] == 1);
    }

    // You may also want to check if the best configuration yields the lowest objective value
    // compared to other random configurations. However, this could be nondeterministic
    // and not suitable for unit tests which should be deterministic. Such a test could
    // be part of a larger integration test suite.
}


// Additional test cases for CalculateObjective function
TEST(QuboSolver, CalculateObjectiveAdditionalCases) {
    double Q[N][N] = {
        {1, -2, 1},
        {-2, 4, -2},
        {1, -2, 1}
    };

    int configuration1[N] = {1, 0, 1};
    double expectedValue1 = 4; // Calculate this by hand for your test case
    double objectiveValue1 = calculateObjective(Q, configuration1);
    EXPECT_DOUBLE_EQ(expectedValue1, objectiveValue1);
    std::cout << "Test Case 1: Configuration {1, 0, 1}, Expected: 4, Actual: " << objectiveValue1 << std::endl;

    int configuration2[N] = {0, 1, 0};
    double expectedValue2 = 2; // Calculate this by hand for your test case
    double objectiveValue2 = calculateObjective(Q, configuration2);
    EXPECT_DOUBLE_EQ(expectedValue2, objectiveValue2);
    std::cout << "Test Case 2: Configuration {0, 1, 0}, Expected: 2, Actual: " << objectiveValue2 << std::endl;
}

// Additional test cases for MonteCarloQUBOSolver function
TEST(QuboSolver, MonteCarloQUBOSolverAdditionalCases) {
    double Q[N][N] = {
        {1, -2, 1},
        {-2, 4, -2},
        {1, -2, 1}
    };
    int numSamples = 1000; // Number of samples for the Monte Carlo simulation
    int bestConfiguration[N];

    monteCarloQUBOSolver(Q, numSamples, bestConfiguration);

    // Ensure the best configuration is a valid binary string
    for (int i = 0; i < N; ++i) {
        EXPECT_TRUE(bestConfiguration[i] == 0 || bestConfiguration[i] == 1);
    }

    // Print the best configuration for visualization
    std::cout << "Best Configuration: ";
    for (int i = 0; i < N; ++i) {
        std::cout << bestConfiguration[i] << " ";
    }
    std::cout << std::endl;
}

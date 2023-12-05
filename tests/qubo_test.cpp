// tests/qubo_test.cpp

#include "gtest/gtest.h"
#include "../src/qubo/qubo.h" // Update the relative path to your quboSolver.h


TEST(QuboSolver, CalculateObjective) {
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
TEST(QuboSolver, MonteCarloQUBOSolver) {
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
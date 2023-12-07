#ifndef QUBO_SOLVER
#define QUBO_SOLVER

namespace QUBO {
const int N = 3; // Example size, you might want this dynamic

double calculateObjective(double** Q, int* configuration, int size);
void monteCarloQUBOSolver(double** Q, int numRoutesPerCar, int numSamples, int* bestConfiguration, int size);
}
#endif
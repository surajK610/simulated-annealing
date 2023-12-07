#ifndef QUBO_SOLVER_H
#define QUBO_SOLVER_H

const int N = 3; // Example size, you might want this dynamic

double calculateObjective(double** Q, int* configuration, int size);
void monteCarloQUBOSolver(double** Q, int numSamples, int* bestConfiguration, int size);
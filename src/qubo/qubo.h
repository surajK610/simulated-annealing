//quboSolver.h

#ifndef QUBO_SOLVER_H
#define QUBO_SOLVER_H

const int N = 3;

double calculateObjective(double Q[N][N], int configuration[N]);

void monteCarloQUBOSolver(double Q[N][N], int numSamples, int bestConfiguration[N]);

#endif

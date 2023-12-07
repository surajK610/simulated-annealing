// main.cpp
#include <cmath>
#include <iostream>

#include "anneal/anneal_omp.hpp"
#include "anneal/anneal.hpp"

#include "schwefel/schwefel.hpp"

void progress(
    void* instance, double cost, float tgen, float tacc, int opt_id, int iter)
{
    printf(
        "bestcost=%1.3e \t tgen=%1.3e \t tacc=%1.3e \t thread=%d\n",
        cost,
        tgen,
        tacc,
        opt_id);
    return ;
}


int main()
{
    srand(0);

    double* x = new double[SCHWEFEL::DIM];
    for (int i = 0; i < SCHWEFEL::DIM; ++i)
        x[i] = drand48();
    double cost = f(nullptr, x);
    printf("Initial cost: %f\n", cost);

    ANNEAL_OMP::SolverOMP solver;
    solver.minimize(SCHWEFEL::DIM, x, f, step, progress, nullptr);

    cost = f(nullptr, x);
    printf("Best cost: %f\nx =\n", cost);

    for (int i = 0; i < SCHWEFEL::DIM; ++i)
        std::cout << 500 * x[i] << " ";
        
    std::cout << std::endl;

    delete[] x;

    return EXIT_SUCCESS;
}

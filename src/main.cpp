// main.cpp
#include <cmath>
#include <iostream>

#include "anneal/anneal_csa.hpp"
#include "anneal/anneal_csa_st.hpp"
#include "anneal/anneal_msa.hpp"
#include "anneal/anneal_msa_st.hpp"
#include "anneal/context.hpp"
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

int main(int argc, char** argv)
{
    int option = 0;
    if (argc > 1)
        option = atoi(argv[1]);

    BaseSolver* solver = nullptr;
    if (option == OPTION_MSA_ST) {
        solver = new MSA_ST::SolverMultipleST();
    } else if (option == OPTION_MSA) {
        solver = new MSA::SolverMultiple();
    } else if (option == OPTION_CSA_ST) {
        solver = new CSA_ST::SolverCoupledST();
    } else if (option == OPTION_CSA) {
        solver = new CSA::SolverCoupled();
    } else {
        std::cout << "Invalid option" << std::endl;
        return EXIT_FAILURE;
    }

    srand(1);

    double* x = new double[SCHWEFEL::DIM];
    for (int i = 0; i < SCHWEFEL::DIM; ++i)
        x[i] = drand48();
    double cost = f(nullptr, x);
    printf("Initial cost: %f\n", cost);
       
    solver->minimize(SCHWEFEL::DIM, x, f, step, progress, nullptr);

    cost = f(nullptr, x);
    printf("Best cost: %f\nx =\n", cost);

    for (int i = 0; i < SCHWEFEL::DIM; ++i)
        std::cout << 500 * x[i] << " ";
        
    std::cout << std::endl;

    delete[] x;

    return EXIT_SUCCESS;
}

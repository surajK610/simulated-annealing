// main.cpp
#include <cmath>
#include <iostream>

#include "anneal/anneal_csa.hpp"
#include "anneal/anneal_csa_st.hpp"
#include "anneal/anneal_msa.hpp"
#include "anneal/anneal_msa_st.hpp"
#include "anneal/context.hpp"
#include "schwefel/schwefel.hpp"
#include "schwefel/schwefel_cuda.hpp"
#include "cuda.h"
#include "curand.h"
#include <sys/time.h>

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
    struct timeval start;
    struct timeval end;

    int optionSA = 0;
    int optionCU = 0;
    if (argc > 1)
        optionSA = atoi(argv[1]);
    if (argc > 2)
        optionCU = atoi(argv[2]);

    BaseSolver* solver = nullptr;
    if (optionSA == OPTION_MSA_ST) {
        solver = new MSA_ST::SolverMultipleST();
    } else if (optionSA == OPTION_MSA) {
        solver = new MSA::SolverMultiple();
    } else if (optionSA == OPTION_CSA_ST) {
        solver = new CSA_ST::SolverCoupledST();
    } else if (optionSA == OPTION_CSA) {
        solver = new CSA::SolverCoupled();
    } else {
        std::cout << "Invalid option for simulated annealing" << std::endl;
        return EXIT_FAILURE;
    }

    gettimeofday(&start, 0);
    if (optionCU == 0) {
        srand(0);
        double* x = new double[SCHWEFEL::DIM];
        for (int i = 0; i < SCHWEFEL::DIM; ++i)
            x[i] = drand48();

        double cost = SCHWEFEL::f(nullptr, x);
        printf("Initial cost: %f\n", cost);
        
        solver->minimize(SCHWEFEL::DIM, x, f, step, progress, nullptr);

        cost = f(nullptr, x);
        printf("Best cost: %f\nx =\n", cost);

        for (int i = 0; i < SCHWEFEL::DIM; ++i)
            std::cout << 500 * x[i] << " ";
            
        std::cout << std::endl;
        delete[] x;
        gettimeofday(&end, 0);
        std::cout << "Took " << (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec) << " microseconds" << std::endl;
        
        return EXIT_SUCCESS;
    } else if (optionCU == 1) {
        srand(0);
        double* x = new double[SCHWEFEL::DIM];
        for (int i = 0; i < SCHWEFEL::DIM; ++i)
            x[i] = drand48();

        dim3 blockSize(256, 1, 1);
        dim3 numBlocks((SCHWEFEL::DIM + blockSize.x - 1) / blockSize.x, 1, 1);

        // setup_kernel<<<numBlocks, blockSize>>>(devStates, time(NULL));

        double cost = f_c(nullptr, x);
        printf("Initial cost: %f\n", cost);

        solver->minimize(SCHWEFEL::DIM, x,  f_c, step_c, progress, nullptr);
        cost = f_c(nullptr, x);

        printf("Best cost: %f\nx =\n", cost);

        for (int i = 0; i <  SCHWEFEL::DIM; ++i)
            std::cout << 500 * x[i] << " ";
            
        std::cout << std::endl;

        delete[] x;
        gettimeofday(&end, 0);
        
        std::cout << "Took " << (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec) << " microseconds" << std::endl;
        return EXIT_SUCCESS;
    } else {
        std::cout << "Invalid option for CUDA" << std::endl;
        return EXIT_FAILURE;
    }
    
}

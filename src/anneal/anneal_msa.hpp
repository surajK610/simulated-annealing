#ifndef ANNEAL_MSA
#define ANNEAL_MSA

#include <cmath>
#include <omp.h>
#include <vector>
#include "context.hpp"

namespace MSA {

class SolverMultiple : public BaseSolver
{
public:
    int m = 4; // number of threads
    int max_iters = 1000000;
    float tgen_initial = 0.01;
    float tgen_schedule = 0.99999;
    float tacc_initial = 0.9;
    float tacc_schedule = 0.01;
    float desired_variance = 0.99;

    SolverMultiple() {  };

    inline int minimize(
        int n, // number of dimensions
        double* x,
        double (*fx)(void*, double*),
        void (*step)(void*, double* y, const double*, float tgen),
        void (*progress)(void*,
                         double cost,
                         float tgen,
                         float tacc,
                         int opt_id,
                         int iter),
        void* instance)
    {
        double fx0 = fx(instance, x);

        SharedStates shared_states(this->m, n, x, fx0);
        float tacc = this->tacc_initial;
        float tgen = this->tgen_initial;
        float gamma = 1;

        omp_lock_t lock;
        omp_init_lock(&lock);

        #pragma omp parallel shared(n, shared_states, tacc, tgen, gamma) num_threads(this->m)
        {
            int k, opt_id = omp_get_thread_num();

            double max_cost = shared_states[0].cost;
            double cost;
            std::vector<double> y(n, double(0));
            float unif, prob;

            #pragma omp for
            for (int iter = 0; iter < this->max_iters; ++iter) {

                step(instance, y.data(), shared_states[opt_id].x.data(), tgen);
                cost = fx(instance, y.data());

                if (cost < shared_states[opt_id].cost) {
                    omp_set_lock(&lock);

                    if (cost < shared_states[opt_id].best_cost) {
                        shared_states[opt_id].best_cost = cost;
                        shared_states[opt_id].best_x = y;
                        if (progress != nullptr)
                            progress(instance, cost, tgen, tacc, opt_id, iter);
                    }

                    shared_states[opt_id].step(y, cost);
                    omp_unset_lock(&lock);
                } else {

                    unif = drand48();
                    prob = std::exp((shared_states[opt_id].cost - max_cost) / tacc) / gamma;
                    if (prob > unif) {
                        omp_set_lock(&lock);
                        shared_states[opt_id].step(y, cost);
                        omp_unset_lock(&lock);
                    }
                }

                tacc += this->tacc_schedule * tacc;
                if (tacc > this->desired_variance)
                    tacc -= this->desired_variance;
                tgen = this->tgen_schedule * tgen;
            }
        }

        int best_ind = 0;
        double best_cost = shared_states[0].best_cost;
        for (int k = 0; k < this->m; ++k) {
            if (shared_states[k].best_cost < best_cost) {
                best_cost = shared_states[k].best_cost;
                best_ind = k;
            }
        }
        Context best_state = shared_states[best_ind];
        for (int i = 0; i < n; ++i)
            x[i] = best_state.best_x[i];

        // Clean up.
        omp_destroy_lock(&lock);

        return 0;
    }
};  // class SolverCoupled
}  // namespace ANNEAL_OMP
#endif


#ifndef ANNEAL
#define ANNEAL

#include <cmath>
#include <vector>
#include <omp.h>

namespace ANNEAL {
class Context
{
public:
    std::vector<double> x;
    std::vector<double> best_x;
    double cost;
    double best_cost;

    /// \param n   The dimension of `x0`.
    /// \param x0  The initial solution guess.
    /// \param fx0 The value of the cost function associated with `x0`.
  
    Context(int n, const double* x0, double fx0)
    {
        this->x = std::vector<double>(x0, x0 + n);
        this->best_x = std::vector<double>(x0, x0 + n);
        this->cost = fx0;
        this->best_cost = fx0;
    }

    /// \param y      The new solution.
    /// \param y_cost The value of the cost function associated with `y`.

    inline void step(std::vector<double> &y, double y_cost)
    {
        this->cost = y_cost;
        this->x.swap(y);
    }
};

class Solver
{
public:
    int m = 1;
    int max_iters = 1000000;
    float tgen_initial = 0.01;
    float tgen_schedule = 0.99999;
    float tacc_initial = 0.9;
    float tacc_schedule = 0.01;
    float desired_variance = 0.99;

    Solver() {  };

    inline int minimize(
        int n,
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

        Context state(n, x, fx0);
        float tacc = this->tacc_initial;
        float tgen = this->tgen_initial;
        float tmp, sum_a, prob_var, gamma = m;

    
        double max_cost = state.cost;
        double cost;
        std::vector<double> y(n, double(0));
        float unif, prob;


        for (int iter = 0; iter < this->max_iters; ++iter) {
            step(instance, y.data(), state.x.data(), tgen);
            cost = fx(instance, y.data());

            if (cost < state.cost) {
                if (cost < state.best_cost) {
                    state.best_cost = cost;
                    state.best_x = y;
                    if (progress != nullptr)
                        progress(instance, cost, tgen, tacc, 0, iter);
                }

                state.step(y, cost);
            } else {
                unif = drand48();
                prob = std::exp((state.cost - max_cost) / tacc) / gamma;
                if (prob > unif) {
                    state.step(y, cost);
                }
            }
            max_cost = state.cost;

            gamma = sum_a = 0.;
    
            tmp = (state.cost - max_cost) / tacc;
            gamma += std::exp(tmp);
            sum_a += std::exp(2.0 * tmp);
            
            prob_var = ((sum_a / (gamma * gamma)) - 1.);

            if (prob_var > this->desired_variance)
                tacc += this->tacc_schedule * tacc;
            else
                tacc -= this->tacc_schedule * tacc;
            tgen = this->tgen_schedule * tgen;
            
        }
        for (int i = 0; i < n; ++i)
          x[i] = state.best_x[i];
        return 0;
    }
};  // class SolverLinear
}  // namespace ANNEAL
#endif


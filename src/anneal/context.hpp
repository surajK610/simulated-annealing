// context.hpp

#ifndef CONTEXT_HPP
#define CONTEXT_HPP

#include <vector>


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

class SharedStates
{
public:
    const int m;
    const int n;

    std::vector<Context> states;

    Context& operator[](int i) {
        return this->states[i];
    }

    const Context& operator[](int i) const {
        return this->states[i];
    }

    /// \param m   The number of threads/shared states.
    /// \param n   The dimension of `x0`.
    /// \param x0  The initial solution guess. Each thread will start from the
    ///            same initial solution.
    /// \param fx0 The value of the cost function associated with `x0`.

    SharedStates(int m, int n, const double* x0, double fx0)
        : m(m), n(n), states(m, Context(n, x0, fx0)) // m different ContextOMP
    {
    }
};  // class SharedStates


enum SolverOption { OPTION_MSA_ST, OPTION_MSA, OPTION_CSA_ST, OPTION_CSA};

class BaseSolver {
public:
    virtual ~BaseSolver() {}
    virtual int minimize(
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
        void* instance) = 0;
};

#endif
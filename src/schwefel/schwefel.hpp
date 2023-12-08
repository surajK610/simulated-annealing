#ifndef SCHWEFEL
#define SCHWEFEL
#include <cmath>
#include <vector>


namespace SCHWEFEL {

const int DIM = 10;

// This defines the Schwefel function. The pointer `instance` can be used to
// pass data to `f`, such as the dimension.
// In this case we don't need to pass anything, so it will be given as NULL.
double f(void* instance, double* x)
{
    double sum = 0.;
    for (int i = 0; i < DIM; ++i)
        sum += 500 * x[i] * std::sin(std::sqrt(std::fabs(500 * x[i])));
    return 418.9829 * DIM - sum;
}

// This function will take a random step from `x` to `y`. The value `tgen`,
// the "generation temperature", determines the variance of the distribution of
// the step. `tgen` will decrease according to fixed a schedule throughout the
// annealing process, which corresponds to a decrease in the variance of steps.
void step(void* instance, double* y, const double* x, float tgen)
{
    int i;
    double tmp;
    for (i = 0; i < DIM; ++i) {
        tmp = std::fmod(x[i] + tgen * std::tan(M_PI * (drand48() - 0.5)), 1.);
        y[i] = tmp;
    }
}
}
#endif
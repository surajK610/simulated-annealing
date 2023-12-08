// schwefel_cuda.hpp
#ifndef SCHWEFEL_CUDA
#define SCHWEFEL_CUDA

#include <vector>
#ifndef __CUDACC__
struct curandState;
#else
#include <curand_kernel.h>
#endif

const int DIM = 10;

// Declaration of CUDA kernels
__global__ void setup_kernel(curandState *state, unsigned long seed);
__global__ void f_kernel(double* x, double* result);
__global__ void step_kernel(double* y, const double* x, float tgen, curandState *state);

// Host functions
void setup_cuda(curandState *dev_states, unsigned long seed);
double f_c(void* instance, double* x);
void step_c(void* instance, double* y, const double* x, float tgen, curandState *dev_states);


#endif // SCHWEFEL_CUDA

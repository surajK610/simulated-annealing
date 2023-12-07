#ifndef SCHWEFEL_CUDA
#define SCHWEFEL_CUDA

#include <cmath>
#include <vector>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

namespace SCHWEFEL_CUDA {

__global__ void f_kernel(double* x, double* result, int dim) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < dim) {
        result[idx] = 500 * x[idx] * sin(sqrt(fabs(500 * x[idx])));
    }
}
__global__ void step_kernel(double* y, const double* x, float tgen, int dim) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < dim) {
        y[idx] = fmod(x[idx] + tgen * tan(M_PI * (drand48() - 0.5)), 1.0);
    }
}

double f_c(void* instance, double* x, int dim) {
    double* dev_x;
    double* dev_result;
    double* host_result = new double[dim];

    cudaMalloc(&dev_x, dim * sizeof(double));
    cudaMalloc(&dev_result, dim * sizeof(double));
    cudaMemcpy(dev_x, x, dim * sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (dim + blockSize - 1) / blockSize;
    f_kernel<<<numBlocks, blockSize>>>(dev_x, dev_result, dim);

    cudaMemcpy(host_result, dev_result, dim * sizeof(double), cudaMemcpyDeviceToHost);

    double sum = 0.;
    for (int i = 0; i < dim; ++i) {
        sum += host_result[i];
    }

    cudaFree(dev_x);
    cudaFree(dev_result);
    delete[] host_result;

    return 418.9829 * dim - sum;
}

void step_c(void* instance, double* y, const double* x, float tgen, int dim) {
    double* dev_x;
    double* dev_y;

    cudaMalloc(&dev_x, dim * sizeof(double));
    cudaMalloc(&dev_y, dim * sizeof(double));
    cudaMemcpy(dev_x, x, dim * sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int numBlocks = (dim + blockSize - 1) / blockSize;
    step_kernel<<<numBlocks, blockSize>>>(dev_y, dev_x, tgen, dim);

    cudaMemcpy(y, dev_y, dim * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(dev_x);
    cudaFree(dev_y);
}
}

#endif
// schwefel_cuda.hpp file

#ifndef SCHWEFEL_CUDA
#define SCHWEFEL_CUDA

#include "schwefel_cuda.hpp"
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <typeinfo>
#include <limits>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>

const int DIM = 10;
const int WARP_SIZE = 32;
const int FULL_MASK = 0xffffffff;
const bool WARP = false;

__global__ void setup_kernel(curandState *state, unsigned long seed) {
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;  
    if (idx < DIM)
      curand_init(seed, idx, 0, &state[idx]); 
}

__global__ void f_kernel(double* x, double* result) {
    size_t idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < DIM) {
      result[idx] = 500 * x[idx] * sin(sqrt(fabs(500 * x[idx])));
    }
}

__global__ void f_kernel_warp(double* x, double* result) {
  // max 1 block
    size_t lane = threadIdx.x % WARP_SIZE;
    int warpid = threadIdx.x/WARP_SIZE;
    int nwarps = blockDim.x/WARP_SIZE;

    double sum = 0.0;
    if (lane < DIM) {
      for (size_t idx = lane + WARP_SIZE*warpid; idx < DIM; idx += WARP_SIZE*nwarps) { // modulus addition
        if (idx < DIM) {
          sum += 500 * x[idx] * sin(sqrt(fabs(500 * x[idx])));
        }
      }
      __syncwarp();
      for (size_t offset = WARP_SIZE/2; offset > 0; offset /= 2) {
        sum += __shfl_down_sync(FULL_MASK, sum, offset);
      }

      __shared__ double s_mem[1024/WARP_SIZE];
      if (lane == 0) {
        s_mem[warpid] = sum;
      }

      __syncthreads(); // sync threads within block
      if (threadIdx.x == 0) { // first lane in first warp
        for (int j = 0; j < nwarps; ++j) {
          result[0] += s_mem[j];
        }   
      }
    }
}

__global__ void step_kernel(double* y, const double* x, float tgen, curandState *state) {
  int idx = threadIdx.x + blockIdx.x * blockDim.x; 
  if (idx < DIM) {
    // curand_init(0, idx, 0, &state[idx]); 
    curandState localState = state[idx];         
    double randVal = curand_uniform(&localState); 
    y[idx] = fmod(x[idx] + tgen * tanf(M_PI * (randVal - 0.5)), 1.0);
    state[idx] = localState;                      
  }
}

double f_c(void* instance, double* x) {
  double *d_x, *d_result;
  double* h_result = new double[DIM];

  cudaMalloc((void**) &d_x, DIM * sizeof(double));
  cudaMalloc((void**) &d_result, DIM * sizeof(double));
  cudaMemcpy(d_x, x, DIM * sizeof(double), cudaMemcpyHostToDevice);

  dim3 blockSize(256, 1, 1);
  dim3 numBlocks((DIM + blockSize.x - 1) / blockSize.x, 1, 1);
  double sum = 0.0;

  if (WARP) {
    if (DIM > 1024) {
      printf("Warp kernel only supports DIM <= 1024\n");
      exit(1);
    }
    
    f_kernel_warp<<<1, 1024>>>(d_x, d_result); // only 1 block, max 1024 threads
    cudaMemcpy(h_result, d_result, 1 * sizeof(double), cudaMemcpyDeviceToHost);
    sum = h_result[0];
  } else {
    f_kernel<<<numBlocks, blockSize>>>(d_x, d_result);
    cudaMemcpy(h_result, d_result, DIM * sizeof(double), cudaMemcpyDeviceToHost);
    for (int i = 0; i < DIM; ++i) {
      sum += h_result[i];
    }
  }
  cudaDeviceSynchronize();
  // fflush(stdout);

  
  delete[] h_result;
  cudaFree(d_x);
  cudaFree(d_result);

  return 418.9829 * DIM - sum;
}

void step_c(void* instance, double* y, const double* x, float tgen) {
  static bool is_setup_done = false;
  static curandState *dev_states;

  double *d_x, *d_y;

  cudaMalloc((void**) &d_x, DIM * sizeof(double));
  cudaMalloc((void**) &d_y, DIM * sizeof(double));
  cudaMemcpy(d_x, x, DIM * sizeof(double), cudaMemcpyHostToDevice);

  dim3 blockSize(256, 1, 1);
  dim3 numBlocks((DIM + blockSize.x - 1) / blockSize.x, 1, 1);

  if (!is_setup_done) {
    cudaMalloc((void**) &dev_states, DIM * sizeof(curandState));
    setup_kernel<<<numBlocks, blockSize>>>(dev_states, time(NULL));
    is_setup_done = true;
    printf("Setup done\n");
  }
  step_kernel<<<numBlocks, blockSize>>>(d_y, d_x, tgen, dev_states);
  cudaDeviceSynchronize();
  cudaMemcpy(y, d_y, DIM * sizeof(double), cudaMemcpyDeviceToHost);

  cudaFree(d_x);
  cudaFree(d_y);

}


#endif
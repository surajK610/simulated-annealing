// schwefel_cuda.cu file

#ifndef SCHWEFEL_CUDA_HPP
#define SCHWEFEL_CUDA_HPP

#include <cmath>
#include <vector>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <iostream>
#include <typeinfo>
#include <limits>

/**
 * Variable Explanations for CUDA Simulated Annealing Algorithm
 *
 * 1. x: Input array for the Schwefel function or the annealing step.
 *    - Type: double*
 *    - Size: 'dim' elements, where each element is a double.
 *    - Role: In the context of the Schwefel function, it represents the point in the function's domain at which the function is evaluated. In the context of the annealing step, it represents the current state of the optimization variable.
 *
 * 2. y: Output array for the annealing step.
 *    - Type: double*
 *    - Size: 'dim' elements, similar to 'x'.
 *    - Role: Stores the new state of the optimization variable after performing the annealing step. Each element of 'y' is computed based on the corresponding element of 'x', with some perturbation determined by 'tgen' and a random number.
 *
 * 3. tgen (generation temperature): Control parameter for the annealing step.
 *    - Type: float
 *    - Role: Determines the variance of the distribution from which the perturbation is drawn during the annealing step. A higher 'tgen' allows for larger jumps in the solution space, facilitating exploration.
 *
 * 4. dim (dimension): The size of the optimization problem.
 *    - Type: int
 *    - Role: Specifies the number of variables in the optimization problem, which directly determines the size of 'x' and 'y' arrays.
 *
 * 5. curandState: Array of states for the CURAND random number generator.
 *    - Type: curandState*
 *    - Size: One state per thread. The total number is usually equal to the total number of threads launched in the kernel.
 *    - Role: These states are used to generate random numbers in a parallel and efficient manner on the GPU. Each thread uses its own state to ensure the independence of random numbers across threads.
 *
 * 6. dev_x, dev_y, dev_result: Device-side counterparts of 'x', 'y', and the intermediate result array.
 *    - Type: double*
 *    - Size: Same as their host-side counterparts ('x', 'y', and an array for storing intermediate results of the Schwefel function).
 *    - Role: These are the device-side (GPU memory) arrays. 'dev_x' and 'dev_y' are used in the annealing step, while 'dev_result' is used to store the output of the Schwefel function computation on the GPU before transferring it back to the host.
 *
 * Note: The size of 'curandState' and the number of threads per block should be carefully chosen to balance efficient use of GPU resources and performance. The block size is often set to a multiple of 32 to align with the warp size of the GPU, enhancing efficiency.
 */

// Setup kernel for initializing CURAND states
__global__ void setup_kernel(curandState *state, unsigned long seed) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;  // Calculate the global thread index
    curand_init(seed, idx, 0, &state[idx]);           // Initialize CURAND state for each thread
}

// Kernel to compute the Schwefel function for each element
__global__ void f_kernel(double* x, double* result, int dim) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;  // Calculate the global thread index
    if (idx < dim) {
        // Compute the Schwefel function for the element at idx
        result[idx] = 500 * x[idx] * sin(sqrt(fabs(500 * x[idx])));
    }
}

// Kernel for the annealing step
__global__ void step_kernel(double* y, const double* x, float tgen, int dim, curandState *state) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;  // Calculate the global thread index
    if (idx < dim) {
        curandState localState = state[idx];          // Get the local CURAND state
        double randVal = curand_uniform(&localState); // Generate a random number
        // Perform the annealing step
        y[idx] = fmod(x[idx] + tgen * tan(M_PI * (randVal - 0.5)), 1.0);
        state[idx] = localState;                      // Update the CURAND state
    }
}

// Host function to compute the Schwefel function
double f_c(void* instance, double* x, int dim, double* dev_x, double* dev_result) {
    // Copy data from host to device
    cudaMemcpy(dev_x, x, dim * sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;  // Block size, chosen for efficient GPU utilization
    int numBlocks = (dim + blockSize - 1) / blockSize;  // Calculate the number of blocks needed
    // Launch the kernel to compute the Schwefel function
    f_kernel<<<numBlocks, blockSize>>>(dev_x, dev_result, dim);

    double* host_result = new double[dim];
    // Copy the result back from device to host
    cudaMemcpy(host_result, dev_result, dim * sizeof(double), cudaMemcpyDeviceToHost);

    double sum = 0.;
    // Accumulate the results to compute the final value of the function
    for (int i = 0; i < dim; ++i) {
        sum += host_result[i];
    }

    delete[] host_result;  // Free the host memory

    return 418.9829 * dim - sum;  // Return the computed value
}

// Host function for the annealing step
void step_c(void* instance, double* y, const double* x, float tgen, int dim, double* dev_x, double* dev_y, curandState *dev_states) {
    // Copy data from host to device
    cudaMemcpy(dev_x, x, dim * sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;  // Block size, chosen for efficient GPU utilization
    int numBlocks = (dim + blockSize - 1) / blockSize;  // Calculate the number of blocks needed
    // Launch the kernel for the annealing step
    step_kernel<<<numBlocks, blockSize>>>(dev_y, dev_x, tgen, dim, dev_states);

    // Copy the result back from device to host
    cudaMemcpy(y, dev_y, dim * sizeof(double), cudaMemcpyDeviceToHost);
}


/////////////////////////////////////////////////////////////////
/// Merged Kernel for Annealing Step and Cost Computation
/////////////////////////////////////////////////////////////////

// Merged Kernel for Annealing Step and Cost Computation
__global__ void step_cost_kernel(double* y, const double* x, float tgen, int dim, curandState* state, double* result) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < dim) {
        // Generate a random number for the annealing step
        curandState localState = state[idx];
        double randVal = curand_uniform(&localState);
        // Perform the annealing step
        double newY = fmod(x[idx] + tgen * tan(M_PI * (randVal - 0.5)), 1.0);
        y[idx] = newY;  // Update the new state
        // Compute the Schwefel function for the new state
        result[idx] = 500 * newY * sin(sqrt(fabs(500 * newY)));
    }
}

// Final computation and reduction of the Schwefel function
__global__ void final_cost_reduction(double* result, int dim) {
    double sum = 0.0;
    for (int i = 0; i < dim; ++i) {
        sum += result[i];
    }
    result[0] = 418.9829 * dim - sum; // Store the final cost in the first element
}

// Host function to perform the annealing step and compute the cost
void step_cost_c(void* instance, double* y, const double* x, float tgen, int dim, double* dev_x, double* dev_y, curandState* dev_states, double* dev_result, cudaStream_t stream) {
    // Copy data from host to device
    cudaMemcpyAsync(dev_x, x, dim * sizeof(double), cudaMemcpyHostToDevice, stream);

    // Calculate the number of blocks and threads
    int blockSize = 256;  // Chosen for efficient GPU utilization
    int numBlocks = (dim + blockSize - 1) / blockSize;

    // Launch the merged kernel for the annealing step and partial cost computation
    step_cost_kernel<<<numBlocks, blockSize, 0, stream>>>(dev_y, dev_x, tgen, dim, dev_states, dev_result);

    // Final reduction to compute the total cost
    final_cost_reduction<<<1, 1, 0, stream>>>(dev_result, dim);

    // Copy the result back from device to host asynchronously
    cudaMemcpyAsync(y, dev_y, dim * sizeof(double), cudaMemcpyDeviceToHost, stream);
    cudaMemcpyAsync((void*)&x[0], dev_result, sizeof(double), cudaMemcpyDeviceToHost, stream); // Copy only the final cost
}

#endif

int main(int argc, char *argv[]) {
    
    // Problem configuration
    int dim = 1024;  // Dimension of the optimization problem
    float tgen = 0.1f;  // Temperature for simulated annealing
    const int iterations = 100;  // Number of iterations for testing

    std::cout << "Data types:\n";
    std::cout << "dim type: " << typeid(dim).name() << "\n";
    std::cout << "tgen type: " << typeid(tgen).name() << "\n";

    // Host memory allocation
    double *h_x = new double[dim];
    double *h_y = new double[dim];
    double f_result;

    // Initialize input array with random values
    std::cout << "Initializing input array...\n";
    // for (int i = 0; i < dim; ++i) {
    //     h_x[i] = static_cast<double>(rand()) / RAND_MAX;
    //     std::cout << "h_x[" << i << "] = " << h_x[i] << "\n";
    // }

    // Device memory allocation
    double *d_x, *d_y, *d_result;
    cudaMalloc(&d_x, dim * sizeof(double));
    cudaMalloc(&d_y, dim * sizeof(double));
    cudaMalloc(&d_result, dim * sizeof(double));

    // CURAND state initialization
    curandState *d_states;
    cudaMalloc(&d_states, dim * sizeof(curandState));

    // Setup CURAND states
    int blockSize = 256;
    int numBlocks = (dim + blockSize - 1) / blockSize;
    setup_kernel<<<numBlocks, blockSize>>>(d_states, time(NULL));
    cudaDeviceSynchronize();  // Wait for kernel completion

    for (int iter = 0; iter < iterations; ++iter) {
        std::cout << "Iteration " << iter + 1 << "...\n";

        // Perform annealing step
        std::cout << "Performing annealing step...\n";
        step_c(nullptr, h_y, h_x, tgen, dim, d_x, d_y, d_states);
        std::cout << "Annealing step completed.\n";

        // Compute Schwefel function
        std::cout << "Computing Schwefel function...\n";
        f_result = f_c(nullptr, h_x, dim, d_x, d_result);
        std::cout << "Schwefel function computed.\n";

        // Print the result
        std::cout << "Result of f_c after iteration " << iter + 1 << ": " << f_result << "\n";
    }

    // Verification (Optional): Check and print some values from h_y
    std::cout << "Verifying output array (h_y):\n";
    for (int i = 0; i < std::min(dim, 10); ++i) {  // Print first 10 values as a sample
        std::cout << "h_y[" << i << "] = " << h_y[i] << "\n";
    }

    // Free resources
    delete[] h_x;
    delete[] h_y;
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_result);
    cudaFree(d_states);

    return 0;
}
#!/bin/bash


echo "Current Working Directory (CWD): $(pwd)"


# Request a GPU partition node and access to 1 GPU
#SBATCH -p 3090-gcondo --gres=gpu:1 --gres-flags=enforce-binding

# Ensures all allocated cores are on the same node
#SBATCH -N 1

# Request 1 CPU core
#SBATCH -n 1

#SBATCH -t 00:05:00
#SBATCH -o outputs/with_gpu.out
#SBATCH -e outputs/with_gpu.err

# Load CUDA module
module load cuda/12.2.2  gcc/10.2   
module load cmake
module load googletest


# Print the files in the current working directory
echo "Files in CWD:"
ls

# Compile CUDA program and run
echo "NVCC Compile:"
nvcc -O2 ./src/schwefel/schwefel_cuda.cu -o ./src/schwefel/schwefel_cuda

echo "Profile:"
chmod +x ./src/schwefel/schwefel_cuda
nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report ./src/schwefel/schwefel_cuda


# ./matrixVectorMult
#!/bin/bash

# SLURM Directives
#SBATCH -p gpu --gres=gpu:1 --gres-flags=enforce-binding
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:05:00
#SBATCH -o outputs/with_gpu.out
#SBATCH -e outputs/with_gpu.err

# Load modules
module load cuda/12.2.2 gcc/10.2
module load cmake
module load googletest

echo "Current Working Directory (CWD): $(pwd)"
echo "Files in CWD:"
ls

# Compile CUDA program
echo "NVCC Compile:"
nvcc -O2 ./src/schwefel/schwefel_cuda.cu -o ./src/schwefel/schwefel_cuda

# Profile the CUDA program
echo "Profile:"
nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report ./src/schwefel/schwefel_cuda

# Uncomment if you want to run your main executable as well
# echo "Run Main Executable:"
# ./$(MAIN_EXECUTABLE)
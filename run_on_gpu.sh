#!/bin/bash

# SLURM Directives
#SBATCH -p gpu --gres=gpu:1 --gres-flags=enforce-binding
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -o outputs/with_gpu.out
#SBATCH -e outputs/with_gpu.err

# Load modules
module load cuda/12.2.2 gcc/10.2
make clean
make
# ./build/main 0 1
# ./build/main 0 1
./build/main 0 1
# ./build/main 3 1
# nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report0 ./build/main 0 1
# nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report1 ./build/main 1 1
# nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report2 ./build/main 2 1
# nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report3 ./build/main 3 1

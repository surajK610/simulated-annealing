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
module load cmake
module load googletest
make clean
make
# ./build/main 0 1
# ./build/main 0 1
# nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report_csa ./build/main 3 1
# nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report_csa ./build/main 3 0

# ./build/main 3 1
nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report_msa_st ./build/main 0 1
nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report_msa ./build/main 1 1
nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report_csa_st ./build/main 2 1
nsys profile --stats=true --force-overwrite=true --output=outputs/gpu_report_csa ./build/main 3 1

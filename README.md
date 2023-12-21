# Simulated Annealing Experiments

## Overview

This repository contains the code for experiments comparing simulated annealing algorithm runtimes with parallelization algorithms. We compared multi-start simulated annealing and coupled simulated annealing on optimization of the Schwefel function (10, 100 dimensional). We experimented with runtime optimization with function inlining, OpenMP, and CUDA kernalization. Please find a more comprehensive overview of our experiments in the Technical Report.

## Reproducing Experiments

We assume a SLURM environment for CUDA experiments.

## Time (Microseconds)
|               | MSA Single Thread | MSA       | CSA Single Thread | CSA      |
|---------------|-------------------|-----------|-------------------|----------|
| CPU           |    3493635        |  876549   |  3884980          | 995846   |
| GPU + CPU     |    1358433        |  344000   |  1563125          | 440082   |
| GPU (Warping) + CPU     |    1208433        |  321389   |  1413214          | 400981   |


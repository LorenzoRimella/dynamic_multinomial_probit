#!/bin/bash
#SBATCH --job-name=filtering_T50_m50
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=g100_usr_pmem
#SBATCH --qos=g100_qos_lprod
#SBATCH --mem=64G
#SBATCH --account=uBG25_EcoRim

# Load necessary modules
module load profile/bioinf
module load autoload r

# Prevent each R worker from spawning additional BLAS/OpenMP threads
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript DMP_cineca/experiments/filtering_likelihood_T50_m50_parallel.R

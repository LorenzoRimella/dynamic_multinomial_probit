#!/bin/bash

#SBATCH --job-name=filtering_likelihood_T10_m10_parallel
#SBATCH --time=72:00:00                   # Maximum wall time (hh:mm:ss)
#SBATCH --nodes=2                         # Request one node
#SBATCH --ntasks=90                        # One task (process) total
#SBATCH --cpus-per-task=1                  # One CPU core per task
#SBATCH --partition=g100_usr_pmem      # Partition (queue) to submit to
#SBATCH --qos=g100_qos_lprod                  # Quality of Service
#SBATCH --mem=64G                          # Memory per node (adjust as needed)
#SBATCH --account=uBG25_EcoRim       # Project account number

# Load necessary modules
module load profile/bioinf
module load autoload r

Rscript DMP_cineca/experiments/filtering_likelihood_T10_m10_parallel.R

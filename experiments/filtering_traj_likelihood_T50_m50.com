#!/bin/bash
#SBATCH --job-name=trajectory_T50_m50
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=g100_usr_pmem
#SBATCH --qos=g100_qos_lprod
#SBATCH --mem=64G
#SBATCH --account=uBG25_EcoRim

# Load necessary modules
module load profile/bioinf
module load autoload r

Rscript DMP_cineca/experiments/filtering_traj_likelihood_T50_m50.R

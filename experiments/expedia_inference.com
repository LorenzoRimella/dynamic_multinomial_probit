#!/bin/bash
#SBATCH --job-name=expedia_EM
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/.."

EXPEDIA_TRAINING_ONLY=1 Rscript expedia_inference.R

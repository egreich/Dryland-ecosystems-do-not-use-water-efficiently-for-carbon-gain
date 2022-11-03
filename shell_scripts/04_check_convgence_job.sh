#!/bin/bash
#SBATCH --job-name=SAM_converge
##SBATCH --workdir
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=100000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu

chmod +x ./scripts/04_check_convergence.R
srun ./scripts/04_check_convergence.R

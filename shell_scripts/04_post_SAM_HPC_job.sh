#!/bin/bash
#SBATCH --job-name=post_SAM
##SBATCH --workdir
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=50000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu

chmod +x ./scripts/04_post_SAM_HPC.R
srun ./scripts/04_post_SAM_HPC.R
chmod +x ./scripts/05_check_convergence.R
srun ./scripts/05_check_convergence.R

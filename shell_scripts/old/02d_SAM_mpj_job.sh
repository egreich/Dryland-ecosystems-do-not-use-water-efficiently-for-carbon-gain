#!/bin/bash
#SBATCH --job-name=SAM_mpj
##SBATCH --workdir
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=150000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu

chmod +x ./scripts/02d_script_SAM_mpj.R
srun ./scripts/02d_script_SAM_mpj.R

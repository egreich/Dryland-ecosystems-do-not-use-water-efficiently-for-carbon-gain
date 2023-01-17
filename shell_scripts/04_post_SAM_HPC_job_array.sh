#!/bin/bash
#SBATCH --job-name=postSAM_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/SAM/log/postSAM_%A_%a.log
#SBATCH --cpus-per-task=1
#SBATCH --time=30:00:00
#SBATCH --mem=100000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=1-40

chmod +x ./scripts/04_post_SAM_HPC.R

site=$(sed -n "$SLURM_ARRAY_TASK_ID"p sitepostEND)
modelv=$(sed -n "$SLURM_ARRAY_TASK_ID"p modelvpostEND)

srun ./shell_scripts/run_post_SAM_job.sh $site $modelv

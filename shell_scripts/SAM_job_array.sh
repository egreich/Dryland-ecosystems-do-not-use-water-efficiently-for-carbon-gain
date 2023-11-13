#!/bin/bash
#SBATCH --job-name=SAM_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/SAM/log/AllSites_11-9_%A_%a.log
##SBATCH --cpus-per-task=3
#SBATCH --time=100:00:00
#SBATCH --mem=40000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
##SBATCH --array=14,15,16
#SBATCH --array=9-13

### %A is monsoon job number %a is interior array index

module load R/4.1.2 # load a specific R version

chmod +x ./scripts/02_script_SAM.R # for permissions
chmod +x ./shell_scripts/run_SAM_job.sh # for permissions

site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedEND)
modelv=$(sed -n "$SLURM_ARRAY_TASK_ID"p modelvEND)

# Run the analysis
srun ./shell_scripts/run_SAM_job.sh $site $seed $modelv

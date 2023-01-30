#!/bin/bash
#SBATCH --job-name=SAM_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/SAM/log/AllSites_%A_%a.log
##SBATCH --cpus-per-task=1
#SBATCH --time=70:00:00
#SBATCH --mem=100000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=49-120
##SBATCH --array=1-144

### %A is monsoon job number %a is interior array index

chmod +x ./scripts/03_script_SAM_HPC.R # for permissions

chain=$(sed -n "$SLURM_ARRAY_TASK_ID"p chainEND)
site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedEND)
modelv=$(sed -n "$SLURM_ARRAY_TASK_ID"p modelvEND)

# Run the analysis
srun ./shell_scripts/run_SAM_job.sh $chain $site $seed $modelv

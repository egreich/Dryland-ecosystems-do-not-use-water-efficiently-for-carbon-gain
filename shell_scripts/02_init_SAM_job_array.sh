#!/bin/bash
#SBATCH --job-name=init_SAM_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/SAM/log/inits_%A_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --time=5:00:00
#SBATCH --mem=40000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=1-8

### %A is monsoon job number %a is interior array index

chmod +x ./shell_scripts/run_init_job.sh # for permissions
chmod +x ./scripts/02_get_inits_SAM_HPC.R # for permissions

site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteinitEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedinitEND)

# Run the analysis
srun ./shell_scripts/run_init_job.sh $site $seed

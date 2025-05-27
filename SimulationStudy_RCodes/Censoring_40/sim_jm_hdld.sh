#!/bin/bash
#SBATCH -p common
#SBATCH --array=1-300
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=16g
#SBATCH --time=24:00:00
#SBATCH --mail-type=all
##SBATCH --mail-user=mohammadsamsul.alam@duke.edu
source activate base
conda activate /hpc/group/dhvi/ma521/miniconda3/envs/env_R
Rscript /cwork/ma521/hdld_surv_jm/RCodes/${SLURM_ARRAY_TASK_ID}.R
conda deactivate



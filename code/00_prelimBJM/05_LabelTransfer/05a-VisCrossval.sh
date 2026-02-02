#!/bin/bash
#SBATCH --job-name=cv
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --output=logs/xval.txt
#SBATCH --error=logs/xval.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=bernie.mulvey@libd.org
#SBATCH --time=2-00:00:00

echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R

## List current modules for reproducibility
module list

## NMF script
Rscript 05a-VisCrossval.R

echo "**** Job ends ****"
date
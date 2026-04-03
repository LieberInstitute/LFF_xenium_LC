#!/bin/bash
#SBATCH --job-name=cv
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --constraint="sapphirerapids"
#SBATCH --output=logs/02_nmfLT63.txt
#SBATCH --error=logs/02_nmfLT63.txt
#SBATCH --mail-type=END
#SBATCH --mail-user=bernie.mulvey@libd.org
#SBATCH --time=4-00:00:00

echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.3.x

## List current modules for reproducibility
module list

## NMF script
Rscript 02-VisToXen_nmfLabelTransfer_k63.R

echo "**** Job ends ****"
date

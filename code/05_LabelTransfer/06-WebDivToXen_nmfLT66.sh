#!/bin/bash
#SBATCH --job-name=lxfr
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --output=logs/06_nmfLT66.txt
#SBATCH --error=logs/06_nmfLT66.txt
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
Rscript 06-WebDivToXen_nmfLT66.R

echo "**** Job ends ****"
date

#!/bin/bash
#SBATCH --mem=60G
#SBATCH --job-name=buildSPE_DAPI
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/buildSPE_DAPI_%a.txt
#SBATCH -e logs/buildSPE_DAPI_%a.txt
#SBATCH --array=2,3,5,11-13,15,17,30%5

set -e


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_CLUSTER_NAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load Modules
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Run code
brnum=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $1}' /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code/image_processing/inputs.txt | awk '{print $1}')
echo "$brnum"

export BRNUM=$brnum
Rscript 01a_build_speNMDAPI.R


## Memeory stat
#sstat -a -o JobID,MaxVMSizeNode,MaxVMSize,AveVMSize,MaxRSS,AveRS S,MaxDiskRead,MaxDiskWrite,AveCPUFreq,TRESUsageInMax -j ${SLURM_JOB_ID}

echo "**** Job ends ****"
date
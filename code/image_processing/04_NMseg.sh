#!/bin/bash
#SBATCH --job-name=NMseg
#SBATCH --mem=60G
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/code/xenium_imageProcessing/logs/NMseg_%a.txt
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/code/xenium_imageProcessing/logs/NMseg_%a.txt
#SBATCH --array=1-29%5

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load MATLAB
module load matlab

matlab -nodesktop -nosplash -r "NMseg"

echo "**** Job ends ****"
date
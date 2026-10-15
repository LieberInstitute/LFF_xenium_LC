#!/bin/bash
#SBATCH --job-name=xLC_appendHE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=125G
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code_newseg/02_QC/logs/writemidresHEDAPIregistSPE.out
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code_newseg/02_QC/logs/writemidresHEDAPIregistSPE.out

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
module load conda_R/4.5

cd /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code_newseg/02_QC/

Rscript --no-save --no-restore 05d-writeSingleSamp_midsizeHE_SPEs.R

echo "**** Job ends ****"
date


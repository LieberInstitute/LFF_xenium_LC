#!/bin/bash
#SBATCH --job-name=BNKS_xLC_mult
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=130G
#SBATCH --array=1-48
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/00_prelim_procdat_BJM/04_Banksy_clusterings/01-banksy_clustruns_outs/logs/ManyBanksy_%a.out
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/00_prelim_procdat_BJM/04_Banksy_clusterings/01-banksy_clustruns_outs/logs/ManyBanksy_%a.err

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.4

## run this script from code dir
cd /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code/00_prelimBJM/04_Banksy_clusterings
Rscript --no-save --no-restore 01b-Many_Banksy.R

echo "**** Job ends ****"
date


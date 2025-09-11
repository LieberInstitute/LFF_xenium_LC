setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/')
library(here)
library(SpatialExperiment)
library(SingleCellExperiment)
library(readxl)
library(tidyverse)

spe_DAPI = readRDS(here("processed-data/xenium/DAPI", "raw_spe.RDS"))
spe_NM_DAPI = readRDS(here("processed-data/xenium/NM_DAPI", "raw_spe.RDS"))

# Step 1: Annotate each SPE with its source
colData(spe_DAPI)$sample_type <- "DAPI"
colData(spe_NM_DAPI)$sample_type <- "NM_DAPI"

# Step 2: Combine the two SPEs
combined_spe <- cbind(spe_DAPI, spe_NM_DAPI)

# Check the result
table(colData(combined_spe)$sample_type)

saveRDS(combined_spe, here("processed-data/xenium/", "raw_combined_spe.RDS"))
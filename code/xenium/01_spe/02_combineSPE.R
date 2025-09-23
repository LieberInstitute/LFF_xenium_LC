setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/')
library(here)
library(SpatialExperiment)
library(SingleCellExperiment)
library(readxl)
library(tidyverse)

#files <- list.files(here('processed-data/xenium/NM_DAPI_rawSPE/'), pattern = "\\.RDS$", full.names = TRUE)
files <- list.files(here('processed-data/xenium/DAPI_rawSPE/'), pattern = "\\.RDS$", full.names = TRUE)

# 1) Read all SPEs and tag brnum
spe_list <- lapply(files, function(f) {
  spe <- readRDS(f)
  br  <- sub(".RDS", " ",basename(f))
  colData(spe)$brnum <- br
  spe
})

combined_spe <- do.call(cbind, spe_list)

#saveRDS(combined_spe, here("processed-data/xenium/NM_DAPI_rawSPE/", "raw_combined_spe.RDS"))
saveRDS(combined_spe, here("processed-data/xenium/DAPI_rawSPE/", "raw_combined_spe.RDS"))
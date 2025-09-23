setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/')
library(here)
library(SpatialExperiment)
library(scran)
library(tidyverse)
library(escheR)
library(scater)
library(scattermore)

DAPI <- readRDS(here('processed-data/xenium/DAPI_rawSPE/raw_combined_spe.RDS'))
NMDAPI <- readRDS('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium/NMDAPI_rawSPE/raw_combined_spe.RDS')

DAPI$sampletype = "DAPI"
NMDAPI$sampletype = "NMDAPI"

spe = cbind(DAPI,NMDAPI)

pdf(here("plots", "xenium1", "02_QC", "DAPIvsNMDAPI.pdf"))
plotColData(spe, y="total_counts", x="sampletype")+geom_scattermore()
plotColData(spe, y="control_probe_counts", x="sampletype")+geom_scattermore()
plotColData(spe, y="unassigned_codeword_counts", x="sampletype")+geom_scattermore()
plotColData(spe, y="cell_area", x="sampletype")+geom_scattermore()
plotColData(spe, y="nucleus_area", x="sampletype")+geom_scattermore()
plotColData(spe, y="transcript_counts", x="sampletype")+geom_scattermore()
dev.off()


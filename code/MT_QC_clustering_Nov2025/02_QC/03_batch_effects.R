
library(tidyverse)
library(escheR)
library(here)
library(SpatialExperiment)
library(scuttle)
library(scattermore)
library(gridExtra)
library(matrixStats)


##########################################################################################
# Run PCA on the Xenium data and colour it by various technical factors
##########################################################################################

# Load the data
spe <- readRDS(here("processed-data", "01_build_spe", "raw_spe_N24.RDS"))
outlier_ids <- read.csv(here("processed-data", "02_xenium_qc", "outlier_ids.csv"))$x

print("propotion of outliers:")
print(length(outlier_ids)/dim(spe)[2])
print("number of unfiltered cells:")
print(dim(spe))
print("number of good quality cells:")
spe <- spe[, -which(colnames(spe) %in% outlier_ids)]
print(dim(spe))

# normalize data
spe <- logNormCounts(spe)

# try pseudobulk
spe_pseudo <- aggregateAcrossCells(spe, ids=spe$BrNum)
spe_pseudo <- logNormCounts(spe_pseudo, size.factors = NULL)

# Run PCA
set.seed(1742)
spe_pseudo <- scater::runPCA(spe_pseudo, ncomponents = 10)


# Get the spe's slide and run date info
slide_id <- unlist(lapply(strsplit(spe_pseudo$Sample, split="/"), "[", 9))
slide_id <- unlist(lapply(strsplit(slide_id, split="__Br"), "[", 1))

spe_pseudo$slide_id <- slide_id

run_date <- unlist(lapply(strsplit(spe_pseudo$Sample, split="/"), "[", 7))
spe_pseudo$run_date <- run_date

pdf(here("plots", "02_xenium_qc", "02_batch_effects_pca", "batch_effects_pca.pdf"), width = 10, height = 10)
scater::plotPCA(spe_pseudo, colour_by = "Dx", ncomponents=4)+
    geom_scattermore()
scater::plotPCA(spe_pseudo, colour_by="PNN", ncomponents=4)
scater::plotPCA(spe_pseudo, colour_by="BrNum", ncomponents=4)
scater::plotPCA(spe_pseudo, colour_by="Sex", ncomponents=4)
scater::plotPCA(spe_pseudo, colour_by="Age", ncomponents=4)
scater::plotPCA(spe_pseudo, colour_by="slide_id", ncomponents=4)
scater::plotPCA(spe_pseudo, colour_by="run_date", ncomponents=4)
dev.off()
 

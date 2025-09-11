setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/')
library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Seurat)
library(scater)
library(cowplot)
library(ggplot2)
library(here)
library(tidyverse)
library(escheR)


##########################################################################################
# Use BANKSY to do integration-free spatial domain analysis as a first pass
##########################################################################################

# Read in the data
spe <- readRDS(here('processed-data/xenium/raw_combined_spe.RDS'))
rowData(spe)$Gene <- rownames(rowData(spe))
colnames(spe) <- paste0(spe$Barcode, "_", spe$sample_id)

empty_cells <- colnames(spe)[colSums(counts(spe)) == 0]
spe <- spe[, colSums(counts(spe)) > 0]

spe <- spe[which(rowData(spe)$Type=="Gene Expression"), ]
print(spe)

#### Banksy parameters ####
lambda <- 0.8
res <- 0.5
compute_agf <- FALSE
use_agf <- FALSE
cnm <- sprintf("clust_M%s_lam%s_k50_res%s", as.numeric(use_agf), lambda, res)


brnums <- unique(spe$sample_type)
if(!file.exists(here("processed-data", "xenium", "03_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))){
  spe <- scuttle::logNormCounts(spe)

  # Create a list of each of the SPEs
  spe_list <- list()

  for (i in 1:length(brnums)) {
    spe_list[[i]] <- spe[, spe$sample_type == brnums[i]]
  }

  k_geom <- 10
  spe_list <- lapply(spe_list, computeBanksy, assay_name = "logcounts", 
                    compute_agf = compute_agf, k_geom = k_geom)

  # Merge the SPEs back 
  spe_joint <- do.call(cbind, spe_list)
  rm(spe_list)
  gc()

  # Run BANKSY PCA

  spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, 
                      lambda = lambda, group = "sample_type", seed = 1000)

  # Run UMAP on the BANKSY PCA embedding
  spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf,  
              lambda = lambda, seed = 1000)

  # BANKSY clustering
  spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)


  clusts <- cbind(colData(spe_joint)[, cnm], rownames(colData(spe_joint)))
  print(head(clusts))
  write.csv(clusts, here("processed-data", "xenium", "03_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))
  
  pdf(here("plots", "xenium", "03_clustering", sprintf("banksy_clustering_%s.pdf", cnm)))
  for (i in 1:length(brnums)){
      sub_spe <- spe_joint[, spe_joint$sample_type == brnums[i]]

      print(head(colData(sub_spe)))

      p <- make_escheR(sub_spe) %>%
          add_ground("Banksy")+
          ggtitle(brnums[[i]])
      print(p)

    }
} else{
  clusts <- read.csv(here("processed-data", "xenium", "03_clustering", sprintf("banksy_clustering_lambda%s_res%s.csv", lambda, res)))
  print(head(clusts))
  colData(spe)[["Banksy"]] <- as.character(clusts$V1)

  # split the spe again and plot each one 
  pdf(here("plots", "xenium", "03_clustering", sprintf("banksy_clustering_%s.pdf", cnm)))
  for (i in 1:length(brnums)){
      sub_spe <- spe[, spe$sample_type == brnums[i]]

      print(head(colData(sub_spe)))

      p <- make_escheR(sub_spe) %>%
          add_ground("Banksy")+
          ggtitle(brnums[[i]])
      print(p)

    }
  dev.off()
}
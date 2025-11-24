setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/')
library(here)
library(SpatialExperiment)
library(scran)
library(tidyverse)
library(escheR)
library(scater)
library(scattermore)

DAPI <- readRDS(here('processed-data/xenium/DAPI_rawSPE/raw_combined_spe.RDS'))
NMDAPI <- readRDS(here('processed-data/xenium/NMDAPI_rawSPE/raw_combined_spe.RDS'))

DAPI$sampletype = "DAPI"
NMDAPI$sampletype = "NMDAPI"

spe = cbind(DAPI,NMDAPI)

pdf(here("plots", "xenium1", "02_QC", "DAPIvsNMDAPI.pdf"))
plotColData(spe, y="total_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
plotColData(spe, y="control_probe_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
plotColData(spe, y="unassigned_codeword_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
plotColData(spe, y="cell_area", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
plotColData(spe, y="nucleus_area", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
plotColData(spe, y="transcript_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
dev.off()

table(spe$sampletype)
#  DAPI NMDAPI 
#507080 548208 

df <- as.data.frame(colData(spe))

library(dplyr)
library(ggplot2)

# 1) Count cells per brnum within each sampletype
df_counts <- df %>%
  mutate(
    brnum = trimws(as.character(brnum)),
    sampletype = trimws(as.character(sampletype))
  ) %>%
  group_by(sampletype, brnum) %>%
  summarise(n_cells = n(), .groups = "drop")

# 2) Box plot by sampletype; points = each brnum (colored)
ggplot(df_counts, aes(x = sampletype, y = n_cells)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(aes(color = brnum), width = 0.15, height = 0, size = 2, alpha = 0.9) +
  labs(x = NULL, y = "Cells per brnum", title = "Cells per brnum by sampletype") +
  scale_y_continuous(labels = scales::comma) +
  theme_classic(base_size = 12)
  
tapply(df$nucleus_area < 4, df$sampletype, sum, na.rm = TRUE)
#DAPI NMDAPI 
#  86  40282

tapply(df$nucleus_area < 5, df$sampletype, sum, na.rm = TRUE)
#DAPI NMDAPI 
#  3158  43315

spe1 <- spe[, spe$nucleus_area > 4]
table(spe1$sampletype)
#  DAPI NMDAPI 
#506994 507926

df1 <- as.data.frame(colData(spe1))

df_counts <- df1 %>%
  mutate(
    brnum = trimws(as.character(brnum)),
    sampletype = trimws(as.character(sampletype))
  ) %>%
  group_by(sampletype, brnum) %>%
  summarise(n_cells = n(), .groups = "drop")

# 2) Box plot by sampletype; points = each brnum (colored)
ggplot(df_counts, aes(x = sampletype, y = n_cells)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(aes(color = brnum), width = 0.15, height = 0, size = 2, alpha = 0.9) +
  labs(x = NULL, y = "Cells per brnum", title = "Cells per brnum by sampletype") +
  scale_y_continuous(labels = scales::comma) +
  theme_classic(base_size = 12)

  pdf(here("plots", "xenium1", "02_QC", "DAPIvsNMDAPI.pdf"))
  plotColData(spe1, y="control_probe_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
  plotColData(spe1, y="control_codeword_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
  plotColData(spe1, y="unassigned_codeword_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
  plotColData(spe1, y="deprecated_codeword_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
  plotColData(spe1, y="cell_area", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
  plotColData(spe1, y="nucleus_area", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
  plotColData(spe1, y="transcript_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
  plotColData(spe1, y="total_counts", x="sampletype", color_by = "nucleus_area")+geom_scattermore()
  dev.off()
  
 library(dplyr)
tapply(df1$control_probe_counts > 5, df1$sampletype, sum, na.rm = TRUE)head
hits <- df1 %>%
  filter(!is.na(control_probe_counts), control_probe_counts > 5) %>%
  select(cell_id, nucleus_area, control_probe_counts, sampletype, brnum)

tapply(df1$unassigned_codeword_counts > 5, df1$sampletype, sum, na.rm = TRUE)
hits <- df1 %>%
  filter(!is.na(unassigned_codeword_counts), unassigned_codeword_counts > 5) %>%
  select(cell_id, nucleus_area, unassigned_codeword_counts, sampletype, brnum)

tapply(df1$control_codeword_counts > 5, df1$sampletype, sum, na.rm = TRUE)
hits <- df1 %>%
  filter(!is.na(control_codeword_counts), control_codeword_counts > 5) %>%
  select(cell_id, nucleus_area, control_codeword_counts, sampletype, brnum)

tapply(df1$deprecated_codeword_counts > 45, df1$sampletype, sum, na.rm = TRUE)
hits <- df1 %>%
  filter(!is.na(deprecated_codeword_counts), deprecated_codeword_counts > 45) %>%
  select(cell_id, nucleus_area, deprecated_codeword_counts, sampletype, brnum)

tapply(df1$transcript_counts > 10000, df1$sampletype, sum, na.rm = TRUE)
hits <- df1 %>%
  filter(!is.na(transcript_counts), transcript_counts > 10000) %>%
  select(cell_id, nucleus_area, transcript_counts, sampletype, brnum)

tapply(df1$total_counts > 10000, df1$sampletype, sum, na.rm = TRUE)
hits <- df1 %>%
  filter(!is.na(total_counts), total_counts > 10000) %>%
  select(cell_id, nucleus_area, total_counts, sampletype, brnum)
  
saveRDS(spe1, here("processed-data/xenium/NMDAPI_rawSPE/", "raw_combined_spe_filter.RDS"))

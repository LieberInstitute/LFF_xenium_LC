setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC')
library(here)

library(S4Vectors)
lc3 <- readRDS(here("processed-data/06-countsOnly_noImgData_QCed_SPE_split_to_tissSections.RDS"))
louvs <- readRDS(here("processed-data/07_featureSelection_dimred_harmony_clustering/04-Initial_louvLeid_res0.5-1-2_clusterings.RDS"))
louvs <- louvs[["HARMONYlmbna_HDG_SVG_2575"]][,c("rn","snnHARMONYlmbna_HDG_SVG_2575_louv_res1")]
colnames(louvs)[2] <- "clusid"
anno <- read.table(here("processed-data/08_validitycheck_25hdg75svg_louv1/10-25hdg75svg_louv1_annots.txt"),header=T)
louvs <- merge(louvs,anno,by="clusid")
louvs <- DataFrame(louvs,row.names=louvs$rn)[colnames(lc3),]
colLabels(lc3) <- louvs$anno 
unique_spd_values <- unique(lc3$label)
lc3_df <- as.data.frame(colData(lc3))
meta_df <- lc3_df[, c("capture_id", "APOE", "Age", "Sex", "Ancestry", "Diagnosis")]
meta_df <- unique(meta_df)

RGB = read.csv(file = here('processed-data/Images/NMvsBGNMvsRGBBGNM/nRGBBGNM.csv'))
NM = read.csv(file = here('processed-data/Images/NMvsBGNMvsRGBBGNM/NMvsBGNM.csv'))
df = merge(RGB, NM, by.x = "fname", by.y = "fname", all.x = TRUE)

df = read.csv(file = here('processed-data/xenium_imageProcessing/NMmetrics.csv'))
# Merge with df by matching df$fname to lc3$base_path
merge_df <- merge(df, meta_df, by.x = "fname", by.y = "capture_id", all.x = TRUE)
merge_df <- merge_df[-c(11, 25, 37, 38), ]
merge_df <- merge_df %>%
  mutate(APOE_carrier = case_when(
    APOE %in% c("E2/ E2", "E2/ E3") ~ "E2",
    APOE %in% c("E3/ E4", "E4/ E4") ~ "E4",
    TRUE ~ NA_character_  # catch-all for other or missing APOE types
  ))
merge_df$slide <- substr(merge_df$fname, 1, nchar(merge_df$fname) - 3)

library(tidyverse)

long_df <- merge_df %>%
  pivot_longer(
    cols = c(RGBNM1, nRGBNM1, lcR1, lcG1, lcB1, BG1, NM1, BG_NM1,
             RGBNM2, nRGBNM2, lcR2, lcG2, lcB2, BG2, NM2, BG_NM2),
    names_to = c(".value", "tissue"),
    names_pattern = "(.*)(1|2)"
  ) %>%
  mutate(tissue = ifelse(tissue == "1", "section1", "section2"))
  
  library(ggplot2)

  group_vars <- c("Sex", "Ancestry", "APOE", "Age", "Diagnosis", "fname", "slide", "APOE_carrier")

  for (var in group_vars) {
  
    # NM vs BG_NM
    p1 <- ggplot(long_df, aes_string(x = "NM", y = "BG_NM", color = var, shape = "tissue")) +
      geom_point(alpha = 0.7) +
      labs(title = paste("NM vs BG_NM colored by", var), x = "NM", y = "BG_NM") +
      theme_bw()
  
    if (var == "Diagnosis") {
      p1 <- p1 + geom_smooth(method = "lm", se = FALSE, aes_string(group = var), linewidth = 0.8)
    }
  
    ggsave(filename = here::here(paste0("plots/NMseg/RGBvsBGNM/NM_vs_BG_NM_by_", var, ".png")),
           plot = p1, width = 6, height = 5)
  
    # NM vs RGBNM
    p2 <- ggplot(long_df, aes_string(x = "NM", y = "RGBNM", color = var, shape = "tissue")) +
      geom_point(alpha = 0.7) +
      labs(title = paste("NM vs RGBNM colored by", var), x = "NM", y = "RGBNM") +
      theme_bw()
  
    if (var == "Diagnosis") {
      p2 <- p2 + geom_smooth(method = "lm", se = FALSE, aes_string(group = var), linewidth = 0.8)
    }
  
    ggsave(filename = here::here(paste0("plots/NMseg/RGBvsBGNM/NM_vs_RGBNM_by_", var, ".png")),
           plot = p2, width = 6, height = 5)
  
    # BG_NM vs RGBNM
    p3 <- ggplot(long_df, aes_string(x = "BG_NM", y = "RGBNM", color = var, shape = "tissue")) +
      geom_point(alpha = 0.7) +
      labs(title = paste("BG_NM vs RGBNM colored by", var), x = "BG_NM", y = "RGBNM") +
      theme_bw()
  
    if (var == "Diagnosis") {
      p3 <- p3 + geom_smooth(method = "lm", se = FALSE, aes_string(group = var), linewidth = 0.8)
    }
  
    ggsave(filename = here::here(paste0("plots/NMseg/RGBvsBGNM/RGBNM_vs_BG_NM_by_", var, ".png")),
           plot = p3, width = 6, height = 5)
  }
  
  
  #### line plots
  
  ggplot(long_df, aes(x = tissue, y = RGBNM, group = fname)) +
    geom_line(alpha = 0.6) +
    geom_point(aes(color = Diagnosis, shape = Sex), size = 2) +
    labs(title = "NM1 to NM2 per sample",
         x = "Section", y = "NM value") +
    theme_bw()
  
	plot_df <- merge_df %>%
	  select(fname, BG_NM1, BG_NM2, NM1, NM2, RGBNM1, RGBNM2) %>%
	  pivot_longer(
	    cols = -fname,
	    names_to = "Measure",
	    values_to = "Value"
	  )

	# Plot
	ggplot(plot_df, aes(x = Measure, y = Value, group = fname, color = fname)) +
	  geom_line(color = "black", linewidth = 0.5) +
	  geom_point(size = 2) +
	  labs(
	    title = "Connecting Dots for NM1, BG_NM1, RGBNM1, NM2, BG_NM2, RGBNM2",
	    x = "Measure", y = "Value"
	  ) +
	  theme_bw() +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1))
	  
	  ## age vs NM
	  ggplot(long_df, aes(x = Age, y = BG_NM)) +
	    geom_point(aes(color = Ancestry, shape = tissue), alpha = 0.7) +
	    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
	    labs(title = "BG_NM vs Age", x = "Age", y = "BG_NM") +
	    theme_bw()
library(SpatialExperiment)
library(SummarizedExperiment)
library(dplyr)
library(readr)

Md <- "/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC"

spe_dir <- file.path(Md, "processed-data/01_spe/NMDAPInewseg_rawSPE")
reg_dir <- file.path(Md, "processed-data/xenium_imageProcessing_new/registrations")

out_dir <- file.path(Md, "processed-data/01_spe/NMDAPInewseg_rawSPE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mpp <- 0.2125

spe_files <- list.files(spe_dir, pattern = "^Br.*\\.RDS$", full.names = TRUE)

for (spe_file in spe_files) {

  brain <- tools::file_path_sans_ext(basename(spe_file))
  message("Processing ", brain)

  spe <- readRDS(spe_file)

  nm_csv <- file.path(reg_dir, brain, "NM_regionprops1.csv")

  if (!file.exists(nm_csv)) {
    warning("Missing NM_regionprops.csv for ", brain)
    next
  }

  nm <- read_csv(nm_csv, show_col_types = FALSE)
  nm <- nm %>%filter(Area >= 20)
  if (nrow(nm) == 0) {
    write_csv(
      tibble(
        brain = character(),
        cell_id = character(),
        NM_object_id = integer(),
        NM_mean_intensity = numeric(),
        NM_max_intensity = numeric(),
        NM_area = numeric(),
        NM_centroid_x_px = numeric(),
        NM_centroid_y_px = numeric(),
        NM_bbox_n_cells = integer(),
        NM_bbox_ambiguous = logical(),
        match_status = character()
      ),
      file.path(out_dir, paste0(brain, "_NM_xenium_cell_matches.csv"))
    )
    next
  }

  coords <- spatialCoords(spe)

  ## SPE coords are microns, convert to registered image pixels
  cell_x_px <- coords[, 1] / mpp
  cell_y_px <- coords[, 2] / mpp
  cell_id <- as.character(colData(spe)$cell_id)

  nm <- nm %>%
    mutate(
      NM_object_id = row_number(),
      #xmin = BBox_X,
      #xmax = BBox_X + BBox_Width,
	  #ymin = BBox_Y,
	  #ymax = BBox_Y + BBox_Height,
	  #xmin = Centroid_X-10,
	  #xmax = Centroid_X+10,
	  #ymin = Centroid_Y-10,
	  #ymax = Centroid_Y+10,
	  xmin = WeightedCentroid_X-10,
	  xmax = WeightedCentroid_X+10,
	  ymin = WeightedCentroid_Y-10,
	  ymax = WeightedCentroid_Y+10
	  
    )

  match_list <- vector("list", nrow(nm))

  for (j in seq_len(nrow(nm))) {

    inside <- which(
      cell_x_px >= nm$xmin[j] &
        cell_x_px <= nm$xmax[j] &
        cell_y_px >= nm$ymin[j] &
        cell_y_px <= nm$ymax[j]
    )

    n_inside <- length(inside)

    if (n_inside == 0) {
      match_list[[j]] <- tibble(
        brain = brain,
        cell_id = NA_character_,
        NM_object_id = nm$NM_object_id[j],
        NM_mean_intensity = nm$MeanIntensity[j],
        NM_max_intensity = nm$MaxIntensity[j],
        NM_min_intensity = nm$MinIntensity[j],
        NM_median_intensity = nm$MedianIntensity[j],
        NM_area = nm$Area[j],
        NM_centroid_x_px = nm$Centroid_X[j],
        NM_centroid_y_px = nm$Centroid_Y[j],
        NM_bbox_n_cells = 0L,
        NM_bbox_ambiguous = FALSE,
        match_status = "no_cell_centroid_in_bbox"
      )
      next
    }

    if (n_inside == 1) {
      matched <- inside
      status <- "unique_match"
    } else {
      d <- sqrt(
        (cell_x_px[inside] - nm$Centroid_X[j])^2 +
          (cell_y_px[inside] - nm$Centroid_Y[j])^2
      )
	  print(d)
      matched <- inside[which.min(d)]
	  matched_dist_px <- min(d)
      status <- "multiple_cells_in_bbox_closest_used"
    }
	
    match_list[[j]] <- tibble(
      brain = brain,
      cell_id = cell_id[matched],
      NM_object_id = nm$NM_object_id[j],
      NM_mean_intensity = nm$MeanIntensity[j],
      NM_max_intensity = nm$MaxIntensity[j],
      NM_min_intensity = nm$MinIntensity[j],
      NM_median_intensity = nm$MedianIntensity[j],
      NM_area = nm$Area[j],
      NM_centroid_x_px = nm$Centroid_X[j],
      NM_centroid_y_px = nm$Centroid_Y[j],
      matched_cell_x_px = cell_x_px[matched],
      matched_cell_y_px = cell_y_px[matched],
      NM_bbox_n_cells = n_inside,
      NM_bbox_ambiguous = n_inside > 1,
      match_status = status,
	  matched_dist_px = matched_dist_px
    )
  }

  match_df <- bind_rows(match_list)
  multi_matches <- match_df %>%
      filter(NM_bbox_n_cells > 1)

  nrow(multi_matches)
  head(multi_matches)
  dim(multi_matches)
  as.data.frame(multi_matches)
  out_csv <- file.path(out_dir, paste0(brain, "_NM_xenium_cell_matches.csv"))
  write_csv(match_df, out_csv)

  message("Saved ", out_csv)
}
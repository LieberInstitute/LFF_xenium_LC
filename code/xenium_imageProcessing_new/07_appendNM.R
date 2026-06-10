library(SpatialExperiment)
library(SummarizedExperiment)
library(dplyr)
library(readr)

Md <- "/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC"

spe_dir <- file.path(Md, "processed-data/01_spe/NMDAPInewseg_rawSPE")
reg_dir <- file.path(Md, "processed-data/xenium_imageProcessing_new/registrations")

out_dir <- file.path(Md, "processed-data/01_spe/NMDAPInewseg_rawSPE_with_NM_intensity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

match_qc_dir <- file.path(out_dir, "NM_match_QC")
dir.create(match_qc_dir, recursive = TRUE, showWarnings = FALSE)

mpp <- 0.2125

spe_files <- list.files(spe_dir, pattern = "^Br.*\\.RDS$", full.names = TRUE)

for (spe_file in spe_files) {

  brain <- tools::file_path_sans_ext(basename(spe_file))
  message("Processing ", brain)

  spe <- readRDS(spe_file)

  nm_csv <- file.path(reg_dir, brain, "NM_regionprops.csv")

  if (!file.exists(nm_csv)) {
    warning("Missing NM_regionprops.csv for ", brain)
    next
  }

  nm <- read_csv(nm_csv, show_col_types = FALSE)

  ## Initialize output columns
  colData(spe)$NM_mean_intensity <- NA_real_
  colData(spe)$NM_max_intensity  <- NA_real_
  colData(spe)$NM_area           <- NA_real_
  colData(spe)$NM_centroid_x_px  <- NA_real_
  colData(spe)$NM_centroid_y_px  <- NA_real_
  colData(spe)$NM_match          <- FALSE
  colData(spe)$NM_bbox_n_cells   <- 0L
  colData(spe)$NM_bbox_ambiguous <- FALSE

  if (nrow(nm) == 0) {
    saveRDS(spe, file.path(out_dir, paste0(brain, ".RDS")))
    next
  }

  ## SPE spatialCoords are in microns, convert to registered image pixels
  coords <- spatialCoords(spe)

  cell_x_px <- coords[, 1] / mpp
  cell_y_px <- coords[, 2] / mpp

  ## Add bbox limits from MATLAB regionprops
  nm <- nm %>%
    mutate(
      NM_object_id = row_number(),
      xmin = BBox_X,
      xmax = BBox_X + BBox_Width,
      ymin = BBox_Y,
      ymax = BBox_Y + BBox_Height,
      n_spe_centroids_in_bbox = NA_integer_,
      matched_cell_index = NA_integer_,
      matched_cell_name = NA_character_,
      match_status = NA_character_
    )

  for (j in seq_len(nrow(nm))) {

    inside <- which(
      cell_x_px >= nm$xmin[j] &
        cell_x_px <= nm$xmax[j] &
        cell_y_px >= nm$ymin[j] &
        cell_y_px <= nm$ymax[j]
    )

    nm$n_spe_centroids_in_bbox[j] <- length(inside)

    if (length(inside) == 0) {
      nm$match_status[j] <- "no_cell_centroid_in_bbox"
      next
    }

    if (length(inside) == 1) {
      matched <- inside
      nm$match_status[j] <- "unique_match"
    }

    if (length(inside) > 1) {
      ## Flag ambiguous bbox, but still assign closest SPE centroid to NM centroid
      d <- sqrt(
        (cell_x_px[inside] - nm$Centroid_X[j])^2 +
          (cell_y_px[inside] - nm$Centroid_Y[j])^2
      )

      matched <- inside[which.min(d)]
      nm$match_status[j] <- "multiple_cells_in_bbox_closest_used"
    }

    nm$matched_cell_index[j] <- matched
    nm$matched_cell_name[j] <- colnames(spe)[matched]

    colData(spe)$NM_mean_intensity[matched] <- nm$MeanIntensity[j]
    colData(spe)$NM_max_intensity[matched]  <- nm$MaxIntensity[j]
    colData(spe)$NM_area[matched]           <- nm$Area[j]
    colData(spe)$NM_centroid_x_px[matched]  <- nm$Centroid_X[j]
    colData(spe)$NM_centroid_y_px[matched]  <- nm$Centroid_Y[j]
    colData(spe)$NM_match[matched]          <- TRUE
    colData(spe)$NM_bbox_n_cells[matched]   <- length(inside)
    colData(spe)$NM_bbox_ambiguous[matched] <- length(inside) > 1
  }

  ## Save updated SPE
  saveRDS(spe, file.path(out_dir, paste0(brain, ".RDS")))

  ## Save per-brain NM matching QC table
  qc_csv <- file.path(match_qc_dir, paste0(brain, "_NM_match_QC.csv"))
  write_csv(nm, qc_csv)

  message("Saved SPE and QC for ", brain)
}
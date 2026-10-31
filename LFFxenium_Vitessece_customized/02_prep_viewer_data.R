## 02_prep_viewer_data.R
## Per-sample data prep for the Vitessce viewer. For one annotated Xenium SpatialExperiment,
## writes to <out_dir>/<sample_id>/:
##   cells.h5ad          AnnData: X = logcounts (cells x genes), layers$counts, obs = cell labels
##                       (categoricals -> Vitessce cell sets), obsm$spatial = centroids in um,
##                       uns$<label>_colors = palette in category order (scanpy convention)
##   domains.geojson     outline polygons for each domain column/level, in the same um frame
##                       (x right, y down, origin = top-left of the Xenium/DAPI image)
##   he_image.tif        the imgData H&E (RGB 8-bit); he_px_per_um in sample_meta.json
##   sample_meta.json    scale factors, image paths/sizes, column roles, polygon parameters
##   qc_domains.png      domain outlines drawn over the H&E and the labelled cells
## and maintains the project-wide <out_dir>/label_palettes.json (seed it with 00_build_palettes.R)
## so a given label gets the same colour in every sample; labels missing from it are appended.
##
## Usage: Rscript 02_prep_viewer_data.R <spe.qs2|.RDS> <dapi_focus.ome.tif> <out_dir>
##   e.g. Rscript 02_prep_viewer_data.R Br0942/Br0942_test_obj.RDS \
##            Br0942/morphology_focus_0000.ome.tif viewer_data

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(data.table)
    library(Matrix)
    library(sf)
    library(anndataR)
    library(jsonlite)
    library(tiff)
    library(xml2)
})

## ---------------------------------------------------------------------------
## Column roles. Patterns are matched against colData names; anything matched by
## domain_pattern is treated as a domain (drawn as outlines), not as a cell label.
## ---------------------------------------------------------------------------
label_cols        <- c("BanksyCluster_Label", "Xfered_snRNA_Label", "Xfered_SRT_Label")
domain_pattern    <- "^Smoothed"
domain_background <- "Other"     # domain level that is not outlined
domain_background_colour <- "#BDBDBD"   # its colour when a domain column is used as a cell set
numeric_obs_cols  <- c("cell_area", "nucleus_area")

## SPE files: qs2 (preferred for the full build) or RDS, by extension
read_spe <- function(p) if (grepl("\\.qs2$", p, ignore.case = TRUE)) qs2::qs_read(p) else readRDS(p)

## Domain-polygon parameters (um unless noted)
tissue_buffer_um  <- 25     # tissue mask = union of discs of this radius around every cell
tissue_gap_um2    <- 20000  # gaps in the tissue mask smaller than this are filled; larger ones
                            # (tears, big vessels) stay as holes that outlines follow
min_island_cells  <- 10     # drop domain pieces containing fewer cells than this
min_hole_cells    <- 10     # fill holes in a domain piece holding fewer cells than this
chaikin_iter      <- 3      # corner-cutting smoothing passes on the Voronoi-edge outline
simplify_tol_um   <- 1.5    # Douglas-Peucker tolerance after smoothing

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------

## Distinct categorical colours (Glasbey set minus white/near-black/near-grey entries),
## used to extend the shared palette file with any labels it has not seen yet.
base_colours <- local({
    g <- unique(c(Polychrome::glasbey.colors(32), Polychrome::palette36.colors(36)))
    rgb <- col2rgb(g)
    lum <- colSums(rgb * c(0.299, 0.587, 0.114)) / 255
    sat <- (apply(rgb, 2, max) - apply(rgb, 2, min)) / 255
    unname(g[lum > 0.15 & lum < 0.9 & sat > 0.25])
})

## Read the shared palette file (seeded with the paper palettes by 00_build_palettes.R), give
## domain columns the outline colour (background level grey), append colours for any label
## levels not yet in the file, write it back, and return named colour vectors per column.
## Cross-process lock (atomic mkdir) so parallel sample runs (build_site.sh) take turns
## updating the shared palette file instead of overwriting each other's additions.
with_file_lock <- function(path, expr, timeout_s = 600) {
    lock <- paste0(path, ".lock")
    t0 <- Sys.time()
    while (!dir.create(lock, showWarnings = FALSE)) {
        if (difftime(Sys.time(), t0, units = "secs") > timeout_s)
            stop("timed out waiting for ", lock, " (remove it if no other run is active)")
        Sys.sleep(0.2)
    }
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
    force(expr)
}

update_palettes <- function(palette_path, obs, label_cols, domain_cols, background) {
    with_file_lock(palette_path, update_palettes_unlocked(palette_path, obs, label_cols, domain_cols, background))
}

update_palettes_unlocked <- function(palette_path, obs, label_cols, domain_cols, background) {
    pal <- if (file.exists(palette_path)) read_json(palette_path, simplifyVector = TRUE) else list()
    if (is.null(pal$domain_outline_colour)) pal$domain_outline_colour <- "#0609C6"
    for (dc in domain_cols) {
        lv <- levels(obs[[dc]])
        pal[[dc]] <- as.list(setNames(ifelse(lv == background, domain_background_colour, pal$domain_outline_colour), lv))
    }
    for (cl in label_cols) {
        have <- if (is.null(pal[[cl]])) character() else unlist(pal[[cl]])
        new_levels <- setdiff(levels(obs[[cl]]), names(have))
        if (length(new_levels)) {
            if (length(have)) warning(cl, ": no palette colour for ", paste(new_levels, collapse = ", "),
                                      "; assigning automatic colours", call. = FALSE)
            free <- setdiff(base_colours, toupper(have))
            if (length(free) < length(new_levels)) free <- c(free, hcl.colors(length(new_levels), "Dynamic"))
            have <- c(have, setNames(free[seq_along(new_levels)], new_levels))
        }
        pal[[cl]] <- as.list(have)
    }
    write_json(pal, palette_path, auto_unbox = TRUE, pretty = TRUE)
    lapply(pal, unlist)
}

## Chaikin corner cutting on a closed ring (matrix with first row == last row)
chaikin_ring <- function(ring, iter) {
    for (k in seq_len(iter)) {
        p <- ring[-nrow(ring), , drop = FALSE]
        q <- p[c(2:nrow(p), 1), , drop = FALSE]
        new <- matrix(t(cbind(0.75 * p + 0.25 * q, 0.25 * p + 0.75 * q)), ncol = 2, byrow = TRUE)
        ring <- rbind(new, new[1, ])
    }
    ring
}

## Flatten any mix of POLYGON / MULTIPOLYGON / GEOMETRYCOLLECTION into an sfc of POLYGONs,
## dropping lines/points that clipping can leave behind.
polys_only <- function(x) {
    split1 <- function(g) {
        if (inherits(g, "POLYGON")) list(g)
        else if (inherits(g, "MULTIPOLYGON")) lapply(unclass(g), st_polygon)
        else if (inherits(g, "GEOMETRYCOLLECTION")) unlist(lapply(g, split1), recursive = FALSE)
        else list()
    }
    st_sfc(unlist(lapply(x, split1), recursive = FALSE))
}

## Drop the holes of a POLYGON for which drop_hole(hole_sfc) is TRUE (vectorised over holes)
fill_holes <- function(poly, drop_hole) {
    rings <- unclass(poly)
    if (length(rings) == 1) return(poly)
    holes <- st_sfc(lapply(rings[-1], function(h) st_polygon(list(h))))
    st_polygon(rings[c(TRUE, !drop_hole(holes))])
}

## Voronoi tessellation of all cells (one polygon per cell, in cell order; unclipped) plus a
## tissue mask used to clip the dissolved domains.
cell_voronoi <- function(xy, buffer_um, gap_um2) {
    pts <- st_geometry(st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y")))
    tissue <- polys_only(st_union(st_buffer(pts, buffer_um, nQuadSegs = 4)))
    tissue <- st_union(st_sfc(lapply(tissue, fill_holes,
                                     drop_hole = function(h) as.numeric(st_area(h)) < gap_um2)))
    vor <- st_collection_extract(st_voronoi(st_union(pts), envelope = st_as_sfc(st_bbox(tissue))), "POLYGON")
    hit <- st_intersects(pts, vor)
    stopifnot(all(lengths(hit) == 1))
    list(cells = vor[unlist(hit)], tissue = tissue, pts = pts)
}

## Outline of the cells in `in_domain`: dissolve their Voronoi polygons, clip to tissue, drop
## small islands, fill holes holding few cells, smooth, simplify, re-clip to tissue.
## Returns an sf of POLYGON parts with the number of domain cells inside each.
domain_outline <- function(vor, in_domain) {
    if (!any(in_domain)) return(NULL)
    dom_pts <- vor$pts[in_domain]
    n_cells_in <- function(g, p) lengths(st_intersects(g, p))
    parts <- polys_only(st_intersection(st_union(vor$cells[in_domain]), vor$tissue))
    parts <- parts[n_cells_in(parts, dom_pts) >= min_island_cells]
    if (length(parts) == 0) return(NULL)
    few_cells <- function(h) n_cells_in(h, vor$pts) < min_hole_cells
    smooth <- lapply(parts, function(poly) {
        st_polygon(lapply(unclass(fill_holes(poly, few_cells)), chaikin_ring, iter = chaikin_iter))
    })
    out <- st_make_valid(st_simplify(st_sfc(smooth), dTolerance = simplify_tol_um))
    out <- polys_only(st_intersection(polys_only(out), vor$tissue))
    n <- n_cells_in(out, dom_pts)
    out <- st_sf(n_cells = n[n >= min_island_cells], geometry = out[n >= min_island_cells])
    out
}

## Donor-level colData (constant within a sample) for viewer labels/titles
donor_cols <- c(brnum = "brnum", apoe = "APOE", sex = "Sex", age = "Age", ancestry = "Ancestry", rin = "Rin")
donor_info <- function(cd) {
    lapply(donor_cols, function(cl) {
        if (!cl %in% names(cd)) return(NULL)
        v <- unique(cd[[cl]])
        if (length(v) != 1) stop("colData column ", cl, " is not constant within the sample")
        v
    })
}

## ---------------------------------------------------------------------------
## Main
## ---------------------------------------------------------------------------
prep_sample <- function(spe_path, dapi_path, out_root) {
    spe <- read_spe(spe_path)
    sid <- unique(spe$sample_id)
    stopifnot(length(sid) == 1)
    out_dir <- file.path(out_root, sid)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("== ", sid, ": ", ncol(spe), " cells x ", nrow(spe), " genes -> ", out_dir)

    xy <- spatialCoords(spe)
    cd <- as.data.table(as.data.frame(colData(spe)), keep.rownames = "cell_id")
    domain_cols <- grep(domain_pattern, names(cd), value = TRUE)
    missing <- setdiff(c(label_cols, numeric_obs_cols), names(cd))
    if (length(missing)) stop("missing colData columns: ", paste(missing, collapse = ", "))

    ## obs: labels as factors (levels sorted, so cell-set order is stable across samples);
    ## domain membership also kept per cell, so it remains usable as a cell set.
    obs <- cd[, c("cell_id", "sample_id", numeric_obs_cols, label_cols, domain_cols), with = FALSE]
    for (cl in c(label_cols, domain_cols)) set(obs, j = cl, value = factor(obs[[cl]], levels = sort(unique(obs[[cl]]))))
    obs[, sample_id := factor(sample_id)]

    ## Shared palettes (labels + domain levels), written back with any new levels appended
    pal <- update_palettes(file.path(out_root, "label_palettes.json"), obs, label_cols, domain_cols, domain_background)

    ## ---- domain polygons -------------------------------------------------
    vor <- cell_voronoi(xy, tissue_buffer_um, tissue_gap_um2)
    dom_list <- list()
    for (dc in domain_cols) for (lv in setdiff(levels(obs[[dc]]), domain_background)) {
        o <- domain_outline(vor, obs[[dc]] == lv)
        if (is.null(o)) next
        o$domain_set <- dc
        o$level <- lv
        o$part <- seq_len(nrow(o))
        o$id <- sprintf("%s|%s|%d", dc, lv, o$part)
        o$area_um2 <- round(as.numeric(st_area(o)), 1)
        o$colour <- pal[[dc]][[lv]]
        dom_list[[length(dom_list) + 1]] <- o
    }
    domains <- do.call(rbind, dom_list)
    domains <- domains[, c("id", "domain_set", "level", "part", "n_cells", "area_um2", "colour")]
    geojson_path <- file.path(out_dir, "domains.geojson")
    if (file.exists(geojson_path)) file.remove(geojson_path)
    st_write(domains, geojson_path, driver = "GeoJSON", quiet = TRUE,
             layer_options = c("COORDINATE_PRECISION=2", "RFC7946=NO"))
    dom_summary <- as.data.table(st_drop_geometry(domains))
    dom_summary[, n_labelled := vapply(seq_len(.N), function(i) sum(obs[[domain_set[i]]] == level[i]), 0L)]
    message("domain outlines:")
    print(dom_summary[, .(parts = .N, cells_inside = sum(n_cells), cells_labelled = n_labelled[1],
                          area_mm2 = round(sum(area_um2) / 1e6, 3)), by = .(domain_set, level)])

    ## ---- AnnData ---------------------------------------------------------
    sce <- SingleCellExperiment(
        assays = list(logcounts = assay(spe, "logcounts"), counts = assay(spe, "counts")),
        colData = DataFrame(as.data.frame(obs[, -"cell_id"]), row.names = obs$cell_id),
        rowData = DataFrame(gene_id = rowData(spe)$ID, row.names = rownames(spe)),
        reducedDims = list(spatial = unname(xy)))   # columns: x, y (um)
    S4Vectors::metadata(sce) <- setNames(lapply(c(label_cols, domain_cols), function(cl) unname(pal[[cl]][levels(obs[[cl]])])),
                                         paste0(c(label_cols, domain_cols), "_colors"))
    ad <- as_AnnData(sce, x_mapping = "logcounts", layers_mapping = c(counts = "counts"),
                     obsm_mapping = c(spatial = "spatial"))
    h5ad_path <- file.path(out_dir, "cells.h5ad")
    if (file.exists(h5ad_path)) file.remove(h5ad_path)
    write_h5ad(ad, h5ad_path, compression = "gzip")

    ## ---- H&E from imgData --------------------------------------------------
    idat <- imgData(spe)
    he_ras <- as.matrix(imgRaster(getImg(spe)))
    he_arr <- array(t(col2rgb(he_ras)) / 255, c(dim(he_ras), 3))
    writeTIFF(he_arr, file.path(out_dir, "he_image.tif"), bits.per.sample = 8L)

    ## ---- DAPI geometry (pixels converted later in Python) ----------------
    ome <- read_xml(readTIFF(dapi_path, payload = FALSE)$description)
    px <- xml_attrs(xml_find_first(ome, "//d1:Pixels", xml_ns(ome)))

    meta <- list(
        sample_id = sid,
        source_spe = normalizePath(spe_path),
        n_cells = ncol(spe), n_genes = nrow(spe),
        donor = donor_info(cd),
        coordinate_units = "micron",
        coordinate_frame = "Xenium: x right, y down, origin at top-left of DAPI image",
        centroid_range_um = list(x = range(xy[, 1]), y = range(xy[, 2])),
        he = list(file = "he_image.tif", image_id = idat$image_id[1],
                  width_px = ncol(he_ras), height_px = nrow(he_ras),
                  px_per_um = idat$scaleFactor[1], um_per_px = 1 / idat$scaleFactor[1]),
        dapi = list(file = normalizePath(dapi_path),
                    width_px = as.integer(px[["SizeX"]]), height_px = as.integer(px[["SizeY"]]),
                    um_per_px = as.numeric(px[["PhysicalSizeX"]])),
        label_columns = label_cols, domain_columns = domain_cols,
        domain_background_level = domain_background, numeric_obs_columns = numeric_obs_cols,
        domain_polygon_params = list(method = "dissolved Voronoi cells of domain-labelled cells",
                                     tissue_buffer_um = tissue_buffer_um, tissue_gap_um2 = tissue_gap_um2, min_island_cells = min_island_cells,
                                     min_hole_cells = min_hole_cells, chaikin_iter = chaikin_iter,
                                     simplify_tol_um = simplify_tol_um),
        palettes = "../label_palettes.json")
    write_json(meta, file.path(out_dir, "sample_meta.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)

    ## ---- QC plot ---------------------------------------------------------
    sf_he <- idat$scaleFactor[1]
    to_he <- function(g) g * sf_he   # um -> H&E pixels (affine scale on sfc)
    n_dc <- length(domain_cols)
    png(file.path(out_dir, "qc_domains.png"), width = 600 * min(n_dc, 3), height = 620 * ceiling(n_dc / 3), res = 100)
    par(mfrow = c(ceiling(n_dc / 3), min(n_dc, 3)), mar = c(0.3, 0.3, 1.8, 0.3))
    lab <- obs$BanksyCluster_Label
    pt_col <- adjustcolor(pal$BanksyCluster_Label[as.character(lab)], 0.7)
    for (dc in domain_cols) {
        plot(as.raster(he_arr * 0.55 + 0.45))
        H <- nrow(he_arr)
        flip <- function(g) { m <- st_geometry(g) * sf_he; m * matrix(c(1, 0, 0, -1), 2) + c(-0.5, H + 0.5) }
        points(xy[, 1] * sf_he - 0.5, H - xy[, 2] * sf_he + 0.5, pch = 16, cex = 0.25, col = pt_col)
        d <- domains[domains$domain_set == dc, ]
        if (nrow(d)) plot(flip(d), add = TRUE, border = pal$domain_outline_colour, lwd = 1.8, col = NA)
        title(sprintf("%s (%d parts)", dc, nrow(d)), cex.main = 0.95)
    }
    invisible(dev.off())

    message("wrote: ", paste(list.files(out_dir), collapse = ", "))
    invisible(list(obs = obs, domains = domains, meta = meta))
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 3) stop("usage: Rscript 02_prep_viewer_data.R <spe.qs2|.RDS> <dapi_focus.ome.tif> <out_dir>")
    prep_sample(args[1], args[2], args[3])
}

## 01_inspect_inputs.R
## Inventory one sample's SpatialExperiment and DAPI image, and check that the coordinate
## systems (cell centroids, imgData H&E, full-res DAPI) line up. Prints the inventory and
## writes three alignment PNGs to <out_dir>.
##
## Usage: Rscript 01_inspect_inputs.R <spe.qs2|.RDS> <morphology_focus_0000.ome.tif> [out_dir]
##   e.g. Rscript 01_inspect_inputs.R Br1039L/Br1039L.qs2 Br1039L/morphology_focus_0000.ome.tif \
##            viewer_data/Br1039L/inspect > viewer_data/Br1039L/inspect/01_inspect_inputs.log

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(data.table)
    library(Matrix)
    library(tiff)
    library(xml2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript 01_inspect_inputs.R <spe.qs2|.RDS> <dapi.ome.tif> [out_dir]")
spe_path  <- args[1]
dapi_path <- args[2]   # Xenium focus image (a single in-focus DAPI plane)
out_dir   <- if (length(args) >= 3) args[3] else dirname(spe_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
read_spe <- function(p) if (grepl("\\.qs2$", p, ignore.case = TRUE)) qs2::qs_read(p) else readRDS(p)

hr <- function(title) cat("\n", strrep("=", 70), "\n== ", title, "\n", strrep("=", 70), "\n", sep = "")

## ---------------------------------------------------------------------------
## 1. SpatialExperiment overview
## ---------------------------------------------------------------------------
spe <- read_spe(spe_path)
hr("SpatialExperiment object")
print(spe)
cat("\nassayNames:", paste(assayNames(spe), collapse = ", "), "\n")
for (a in assayNames(spe)) cat(sprintf("  %-10s class = %s\n", a, class(assay(spe, a))[1]))
cat("reducedDimNames:", paste(reducedDimNames(spe), collapse = ", "), "\n")
for (r in reducedDimNames(spe)) cat(sprintf("  %-12s %d x %d\n", r, nrow(reducedDim(spe, r)), ncol(reducedDim(spe, r))))
cat("altExpNames:", paste(altExpNames(spe), collapse = ", "), "\n")
for (a in altExpNames(spe)) { cat("  altExp", a, ":\n"); print(altExp(spe, a)) }
cat("sample_id(s):", paste(unique(spe$sample_id), collapse = ", "), "\n")

## ---------------------------------------------------------------------------
## 2. colData columns: type, # unique values, NA count, example values
## ---------------------------------------------------------------------------
hr("colData columns")
cd <- as.data.table(as.data.frame(colData(spe)), keep.rownames = "cell_id")
cd_summary <- rbindlist(lapply(names(cd), function(nm) {
    v <- cd[[nm]]
    nu <- uniqueN(v)
    ex <- if (is.numeric(v) && nu > 20) {
        sprintf("range [%s, %s]", signif(min(v, na.rm = TRUE), 5), signif(max(v, na.rm = TRUE), 5))
    } else {
        paste(head(sort(unique(as.character(v))), 6), collapse = ", ")
    }
    data.table(column = nm, class = class(v)[1], n_unique = nu, n_NA = sum(is.na(v)),
               examples = substr(ex, 1, 90))
}))
print(cd_summary, nrows = 500, class = FALSE)

## Full level tables for the low-cardinality (label-like) columns
hr("Label-like columns (<= 60 levels): counts per level")
label_cols <- cd_summary[n_unique <= 60 & n_unique > 1 & class %in% c("factor", "character", "integer", "logical"), column]
for (nm in label_cols) {
    cat("\n--", nm, "--\n")
    print(cd[, .N, by = nm][order(-N)], nrows = 100, class = FALSE)
}

## ---------------------------------------------------------------------------
## 3. rowData + logcounts sanity
## ---------------------------------------------------------------------------
hr("rowData / features")
rd <- as.data.table(as.data.frame(rowData(spe)), keep.rownames = "feature")
print(head(rd, 5), class = FALSE)
cat("n features:", nrow(spe), "\n")
if ("logcounts" %in% assayNames(spe)) {
    lc <- assay(spe, "logcounts")
    cat("logcounts class:", class(lc)[1], " dim:", dim(lc),
        " density:", signif(Matrix::nnzero(lc) / prod(dim(lc)), 3),
        " value range:", signif(range(lc@x), 4), "\n")
    top <- sort(Matrix::rowMeans(lc), decreasing = TRUE)
    cat("top-10 genes by mean logcounts:\n"); print(signif(head(top, 10), 3))
}

## ---------------------------------------------------------------------------
## 4. spatialCoords
## ---------------------------------------------------------------------------
hr("spatialCoords")
xy <- spatialCoords(spe)
cat("spatialCoordsNames:", paste(spatialCoordsNames(spe), collapse = ", "), "\n")
xy_dt <- data.table(axis = colnames(xy),
                    min = apply(xy, 2, min), max = apply(xy, 2, max),
                    span = apply(xy, 2, function(v) diff(range(v))))
print(xy_dt, class = FALSE)

## ---------------------------------------------------------------------------
## 5. imgData (embedded image + scale factors)
## ---------------------------------------------------------------------------
hr("imgData")
idat <- imgData(spe)
print(idat)
for (i in seq_len(nrow(idat))) {
    img <- getImg(spe, sample_id = idat$sample_id[i], image_id = idat$image_id[i])
    ras <- imgRaster(img)
    sf  <- idat$scaleFactor[i]
    cat(sprintf("\nimage '%s' (%s): class %s, %d x %d px (W x H), scaleFactor = %s\n",
                idat$image_id[i], idat$sample_id[i], class(img)[1], ncol(ras), nrow(ras), sf))
    cat(sprintf("  full-res extent implied by scaleFactor: %.1f x %.1f (coord units)\n",
                ncol(ras) / sf, nrow(ras) / sf))
    cat(sprintf("  cell centroids x scaleFactor -> pixel range x [%.1f, %.1f], y [%.1f, %.1f]\n",
                min(xy[, 1]) * sf, max(xy[, 1]) * sf, min(xy[, 2]) * sf, max(xy[, 2]) * sf))
}

## (Br0942.tiff was a mis-saved export of the imgData image and is not used.)
ras <- as.matrix(imgRaster(getImg(spe)))
emb <- array(t(col2rgb(ras)) / 255, c(dim(ras), 3))   # H x W x RGB, 0-1

## ---------------------------------------------------------------------------
## 7. DAPI OME-TIFF: OME-XML + pyramid levels
## ---------------------------------------------------------------------------
hr(paste0(dapi_path, " (DAPI)"))
## tiff::readTIFF(payload = FALSE) reads only tags (no pixel decoding), so this is cheap.
dapi_info <- readTIFF(dapi_path, payload = FALSE)
ome <- read_xml(dapi_info$description)
ns  <- xml_ns_rename(xml_ns(ome), d1 = "ome")
px  <- xml_find_first(ome, "//ome:Pixels", ns)
pixel_attrs <- xml_attrs(px)
print(as.data.table(as.list(pixel_attrs)), class = FALSE)
chan <- xml_attrs(xml_find_all(ome, "//ome:Channel", ns))
cat("channels:", paste(vapply(chan, `[[`, "", "Name"), collapse = ", "), "\n")
cat("instrument:", paste(xml_attrs(xml_find_first(ome, "//ome:Microscope", ns)), collapse = " "), "\n")
um_per_px <- as.numeric(pixel_attrs[["PhysicalSizeX"]])

## Pyramid levels in Xenium OME-TIFFs are stored as SubIFDs (tag 330) of each plane,
## which the tiff package does not expose (and RBioFormats is not installed), so use a
## minimal TIFF/BigTIFF IFD reader. Offsets exceed 2^31, so unsigned ints are assembled
## from raw bytes as doubles.
read_uint <- function(con, nbytes, n = 1) {
    b <- matrix(as.numeric(readBin(con, "raw", nbytes * n)), nrow = nbytes)
    colSums(b * 256^(0:(nbytes - 1)))
}
read_ifd <- function(con, offset, big, keep = c(256, 257, 259, 322, 323, 324, 325, 330)) {
    seek(con, offset)
    n <- read_uint(con, if (big) 8 else 2)
    inline_bytes <- if (big) 8 else 4
    type_size <- c(`1` = 1, `3` = 2, `4` = 4, `16` = 8, `18` = 8)
    tags <- list()
    for (k in seq_len(n)) {
        tag   <- read_uint(con, 2)
        type  <- read_uint(con, 2)
        count <- read_uint(con, if (big) 8 else 4)
        here  <- seek(con)
        sz    <- type_size[as.character(type)]
        if (tag %in% keep && !is.na(sz)) {
            if (count * sz > inline_bytes) seek(con, read_uint(con, inline_bytes))
            tags[[as.character(tag)]] <- read_uint(con, sz, count)
        }
        seek(con, here + inline_bytes)
    }
    attr(tags, "next") <- read_uint(con, inline_bytes)   # offset of next IFD in the chain (0 = end)
    tags
}
con <- file(dapi_path, "rb")
hdr <- read_uint(con, 2, 2)
big <- hdr[2] == 43
if (big) invisible(read_uint(con, 2, 2))   # bytesize + reserved
first_ifd <- read_uint(con, if (big) 8 else 4)
## main IFD chain = one IFD per z-plane; each plane carries its own pyramid as SubIFDs
planes <- list(); off <- first_ifd
while (off > 0) { planes[[length(planes) + 1]] <- read_ifd(con, off, big); off <- attr(planes[[length(planes)]], "next") }
## pyr[[z]][[level + 1]] = IFD for plane z at that pyramid level
pyr <- lapply(planes, function(p) c(list(p), lapply(p[["330"]], function(o) read_ifd(con, o, big))))
close(con)
ifd0 <- planes[[1]]
sub_ifds <- pyr[[1]]
levels <- rbindlist(lapply(seq_along(sub_ifds), function(j) {
    s <- sub_ifds[[j]]
    data.table(level = j - 1L, width = s[["256"]], height = s[["257"]],
               n_tiles = length(s[["324"]]), compression = s[["259"]])
}))
levels[, downsample := round(ifd0[["256"]] / width, 2)]
levels[, um_per_px := um_per_px * ifd0[["256"]] / width]
cat(sprintf("\n%s, %d z-planes in OME-XML, %d plane IFDs found, tile %d x %d, compression %d (34712 = JPEG 2000)\n",
            if (big) "BigTIFF" else "TIFF", as.integer(pixel_attrs[["SizeZ"]]), length(planes),
            ifd0[["322"]], ifd0[["323"]], ifd0[["259"]]))
cat("pyramid (SubIFDs; identical for every z-plane):\n")
print(levels, class = FALSE)
cat(sprintf("full-res DAPI extent: %.1f x %.1f um\n",
            ifd0[["256"]] * um_per_px, ifd0[["257"]] * um_per_px))

## ---------------------------------------------------------------------------
## 8. Coordinate-system cross-check
## ---------------------------------------------------------------------------
hr("Coordinate cross-check")
sf_lowres <- idat$scaleFactor[1]
cat(sprintf("cell centroid extent: x [%.1f, %.1f], y [%.1f, %.1f]\n",
            min(xy[, 1]), max(xy[, 1]), min(xy[, 2]), max(xy[, 2])))
cat(sprintf("DAPI full-res extent (um): %.1f x %.1f\n",
            ifd0[["256"]] * um_per_px, ifd0[["257"]] * um_per_px))
cat(sprintf("lowres H&E extent / scaleFactor:  %.1f x %.1f\n",
            ncol(ras) / sf_lowres, nrow(ras) / sf_lowres))
cat(sprintf("=> centroids in um map to DAPI px x [%.0f, %.0f], y [%.0f, %.0f] (image %.0f x %.0f)\n",
            min(xy[, 1]) / um_per_px, max(xy[, 1]) / um_per_px,
            min(xy[, 2]) / um_per_px, max(xy[, 2]) / um_per_px,
            ifd0[["256"]], ifd0[["257"]]))
cat(sprintf("=> lowres scaleFactor %.6f px/um  = %.4f um/px = DAPI um/px x %.3f\n",
            sf_lowres, 1 / sf_lowres, (1 / sf_lowres) / um_per_px))

## ---------------------------------------------------------------------------
## 9. Visual check: centroids over the imgData H&E and a DAPI pyramid level
## ---------------------------------------------------------------------------
hr("Overlay plots")
## Decode one pyramid level: each tile is a bare JPEG 2000 codestream; write it out and
## let GDAL's OpenJPEG driver (via sf) decode it, then stitch tiles.
decode_level <- function(path, ifd, tmpdir = tempdir()) {
    tw <- ifd[["322"]]; th <- ifd[["323"]]; W <- ifd[["256"]]; H <- ifd[["257"]]
    ntx <- ceiling(W / tw); nty <- ceiling(H / th)
    out <- matrix(0, nrow = nty * th, ncol = ntx * tw)
    con <- file(path, "rb"); on.exit(close(con))
    for (k in seq_along(ifd[["324"]])) {
        seek(con, ifd[["324"]][k])
        j2k <- file.path(tmpdir, sprintf("tile%03d.j2k", k)); tif <- sub("j2k$", "tif", j2k)
        writeBin(readBin(con, "raw", ifd[["325"]][k]), j2k)
        sf::gdal_utils("translate", j2k, tif, options = c("-of", "GTiff", "-q"))
        tile <- readTIFF(tif, as.is = TRUE)
        r0 <- ((k - 1) %/% ntx) * th; c0 <- ((k - 1) %% ntx) * tw
        out[r0 + seq_len(nrow(tile)), c0 + seq_len(ncol(tile))] <- tile
    }
    out[seq_len(H), seq_len(W)]
}
## If a z-stack is supplied (e.g. morphology.ome.tif), take a max-intensity projection across
## z; for the focus image there is a single plane and this is a no-op.
dapi_level <- levels[n_tiles <= 4][1, level]           # ~1.4k px wide, 4 tiles
per_z <- lapply(pyr, function(p) decode_level(dapi_path, p[[dapi_level + 1]]))
if (length(per_z) > 1) cat("mean DAPI intensity per z-plane:", round(vapply(per_z, mean, 0), 1), "\n")
dapi_img <- Reduce(pmax, per_z)
dapi_umpp <- levels[level == dapi_level, um_per_px]
cat(sprintf("decoded DAPI level %d: %d x %d px, %.3f um/px, intensity range %s\n",
            dapi_level, ncol(dapi_img), nrow(dapi_img), dapi_umpp, paste(range(dapi_img), collapse = "-")))
q <- quantile(dapi_img[dapi_img > 0], c(0.01, 0.995))
dapi_disp <- pmin(pmax((dapi_img - q[1]) / (q[2] - q[1]), 0), 1)

lab <- factor(spe$BanksyCluster_Label)
pal <- setNames(hcl.colors(nlevels(lab), "Dark 3"), levels(lab))
pt_col <- adjustcolor(pal[lab], alpha.f = 0.6)

overlay <- function(img, px_x, px_y, main, cex = 0.35) {
    plot(as.raster(img)); title(main, cex.main = 1)
    ## plot.raster puts pixel (col c, row r) at x = c - 0.5, y = H - r + 0.5
    H <- if (is.matrix(img)) nrow(img) else dim(img)[1]
    points(px_x - 0.5, H - px_y + 0.5, pch = 16, cex = cex, col = pt_col)
}
png(file.path(out_dir, "01_alignment_check.png"), width = 1700, height = 880, res = 110)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
overlay(emb, xy[, 1] * sf_lowres, xy[, 2] * sf_lowres,
        sprintf("imgData lowres H&E  (px = um x %.4f)", sf_lowres))
overlay(dapi_disp, xy[, 1] / dapi_umpp, xy[, 2] / dapi_umpp,
        sprintf("DAPI, level %d  (px = um / %.2f)", dapi_level, dapi_umpp))
invisible(dev.off())

## zoomed crop to judge cell-level registration on DAPI (level 2, ~0.85 um/px)
zl <- 2L
z_umpp <- levels[level == zl, um_per_px]
z_ifd <- sub_ifds[[zl + 1]]
## pick the tile containing the median cell
cx <- median(xy[, 1]) / z_umpp; cy <- median(xy[, 2]) / z_umpp
ntx <- ceiling(z_ifd[["256"]] / z_ifd[["322"]])
k <- (cy %/% z_ifd[["323"]]) * ntx + (cx %/% z_ifd[["322"]]) + 1
tile <- Reduce(pmax, lapply(pyr, function(p) {
    one <- p[[zl + 1]]; one[["324"]] <- one[["324"]][k]; one[["325"]] <- one[["325"]][k]
    one[["256"]] <- one[["322"]]; one[["257"]] <- one[["323"]]
    decode_level(dapi_path, one)
}))
c0 <- ((k - 1) %% ntx) * z_ifd[["322"]]; r0 <- ((k - 1) %/% ntx) * z_ifd[["323"]]
q <- quantile(tile[tile > 0], c(0.01, 0.995)); tile <- pmin(pmax((tile - q[1]) / (q[2] - q[1]), 0), 1)
png(file.path(out_dir, "01_alignment_zoom_dapi.png"), width = 1100, height = 1150, res = 110)
par(mar = c(0.5, 0.5, 2, 0.5))
overlay(tile, xy[, 1] / z_umpp - c0, xy[, 2] / z_umpp - r0, cex = 0.8,
        sprintf("DAPI, level %d tile %d (%.0f um square) + centroids", zl, k, z_ifd[["322"]] * z_umpp))
invisible(dev.off())
cat("wrote 01_alignment_check.png and 01_alignment_zoom_dapi.png\n")

## ---------------------------------------------------------------------------
## 10. Image-to-image check: imgData H&E vs DAPI, both mapped through um
## ---------------------------------------------------------------------------
hr("H&E (imgData) vs DAPI registration")
## Nearest-neighbour resample the DAPI level onto the H&E grid: H&E pixel centre -> um -> DAPI px
H <- nrow(emb); W <- ncol(emb)
rr <- pmin(pmax(floor(((seq_len(H) - 0.5) / sf_lowres) / dapi_umpp) + 1, 1), nrow(dapi_img))
cc <- pmin(pmax(floor(((seq_len(W) - 0.5) / sf_lowres) / dapi_umpp) + 1, 1), ncol(dapi_img))
dapi_on_he <- dapi_img[rr, cc]
## Tissue masks: smoothed DAPI signal vs non-background H&E; residual shift = argmax IoU
box <- function(m, k = 5) {
    f <- rep(1 / k, k)
    apply(t(apply(m, 1, stats::filter, f, circular = TRUE)), 2, stats::filter, f, circular = TRUE)
}
dsm <- box(dapi_on_he)
dapi_mask <- dsm > quantile(dsm, 0.35)
lum <- (emb[, , 1] + emb[, , 2] + emb[, , 3]) / 3
he_mask <- lum < 0.85 & lum > 0.05
iou <- function(dy, dx) {
    a <- he_mask[max(1, 1 + dy):min(H, H + dy), max(1, 1 + dx):min(W, W + dx)]
    b <- dapi_mask[max(1, 1 - dy):min(H, H - dy), max(1, 1 - dx):min(W, W - dx)]
    sum(a & b) / sum(a | b)
}
shifts <- CJ(dy = -10:10, dx = -10:10)
shifts[, IoU := mapply(iou, dy, dx)]
cat(sprintf("tissue-mask IoU at zero shift: %.3f; best shift (H&E px, 1 px = %.2f um):\n",
            shifts[dy == 0 & dx == 0, IoU], 1 / sf_lowres))
print(shifts[order(-IoU)][1:3], class = FALSE)

q <- quantile(dapi_on_he, c(0.5, 0.995))
dn <- pmin(pmax((dapi_on_he - q[1]) / (q[2] - q[1]), 0), 1)
comp <- emb * 0.6
comp[, , 2] <- pmax(comp[, , 2], dn)
png(file.path(out_dir, "01_alignment_he_vs_dapi.png"), width = 1600, height = 820, res = 110)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2, 0.5))
plot(as.raster(comp)); title("imgData H&E (dimmed) + DAPI (green)")
z <- 300:520
plot(as.raster(comp[z, z, ])); title(sprintf("zoom: H&E px %d-%d", min(z), max(z)))
invisible(dev.off())
cat("wrote 01_alignment_he_vs_dapi.png\n")

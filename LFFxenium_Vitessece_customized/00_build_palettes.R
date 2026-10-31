## 00_build_palettes.R
## Seed the project-wide viewer palette file (viewer_data/label_palettes.json) from the paper
## palettes, renaming their display names (with ggplot markdown) to the label values in the SPEs.
## 02_prep_viewer_data.R keeps these colours and only appends colours for labels not listed here
## (e.g. Xfered_snRNA_Label, which has no paper palette).
##
## Usage: Rscript 00_build_palettes.R [out_dir = viewer_data]

suppressPackageStartupMessages(library(jsonlite))

out_dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "viewer_data"
dir.create(out_dir, showWarnings = FALSE)

## paper palette name -> SPE label value
srt_names <- c(
    "LC"                                  = "LC",
    "*ACHE*, *SERT*, LC, Oligo Mixed"     = "ACHE_SERT_LC_WM",
    "*ACHE*, *SERT*, Neuropeptide Mixed"  = "ACHE_SERT_Npep",
    "Astro"                               = "Astro",
    "Astro-Oligo Mixed"                   = "Astro_Oligo",
    "Vascular"                            = "Vascular",
    "Oligo"                               = "Oligo",
    "Oligo-Astro Mixed"                   = "Oligo_Astro",
    "Oligo and Other Types Mixed"         = "Oligo_Complex")

banksy_names <- c(
    "Blood"            = "Blood",
    "OPC"              = "OPC",
    "Mixed Npept."     = "Mixed Neuropept.",
    "Astro *FGFR3*"    = "Astro_FGFR3",
    "Microglia"        = "Microglia",
    "NK/T Cells"       = "NK_and_T_cells",
    "Oligodend."       = "Oligodendrocytes",
    "GABA *5HTR2A*"    = "Inhib_5HTR2A",
    "Astro *RYR3*"     = "Astro_RYR3",
    "Astro *SERPINA3*" = "Astro_SERPINA3",
    "Ependyma"         = "Ependyma",
    "Vasc. SMC"        = "Vasc. SMC",
    "LCNE"             = "LCNE",
    "Endothel."        = "Endothelium",
    "5HTergic"         = "Serotonergic")
banksy_extra <- c(DROPPED_CLUSTER = "#4D4D4D")   # not in the paper palette

rename_pal <- function(pal, name_map) {
    unknown <- setdiff(names(pal), names(name_map))
    if (length(unknown)) stop("palette entries without a label mapping: ", paste(unknown, collapse = "; "))
    setNames(unname(pal), name_map[names(pal)])
}

srt    <- rename_pal(readRDS("finalpal_somuchfornorainbows.RDS"), srt_names)
banksy <- c(rename_pal(readRDS("banksypal.RDS"), banksy_names), banksy_extra)

pal <- list(
    Xfered_SRT_Label    = as.list(srt),
    BanksyCluster_Label = as.list(banksy),
    ## all domain outlines are drawn in the SRT palette's LC blue
    domain_outline_colour = unname(srt[["LC"]]))

palette_path <- file.path(out_dir, "label_palettes.json")
write_json(pal, palette_path, auto_unbox = TRUE, pretty = TRUE)
message("wrote ", palette_path)
str(pal)

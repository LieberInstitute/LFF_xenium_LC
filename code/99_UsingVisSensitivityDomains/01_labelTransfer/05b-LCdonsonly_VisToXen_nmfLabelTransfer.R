setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC')

# Load libraries
library(data.table)
library(Matrix) # otherwise rcppml gives error 'no item called "package:Matrix" on the search list'
library(RcppML)
library(SingleCellExperiment)
library(SpatialExperiment)
library(nmfLabelTransfer)

lcv <- readRDS("../LFF_spatial_LC/processed-data/06-LCposDonorsOnly_countsOnly_noImgData_QCed_SPE_split_to_tissSections.RDS")

louvs <- fread("../LFF_spatial_LC/processed-data/16_sensitivityanalysis_reproc_wo3donors/Sens1_LCsampsonly_featureSelection_dimred_harmony_clustering/04_2-LConly_HDG_SVG_2575_Louv1_clustering.txt")
setnames(louvs,2,"clusid")

finalann <- fread("../LFF_spatial_LC/processed-data/16_sensitivityanalysis_reproc_wo3donors/Sens1_LCsampsonly_featureSelection_dimred_harmony_clustering/08_2-25hdg75svg_louv1_annots.txt")

louvs <- merge.data.table(louvs,finalann,by="clusid")
louvs <- DataFrame(louvs,row.names=louvs$rn)[colnames(lcv),]

ncol(lcv) # sanity check-120k instead for reprocessed data instead of 130somethingk
colLabels(lcv) <- louvs$anno

rm(vdom,vanno)

lcv$brnum[lcv$brnum=="Br6119(re-dis)"] <- "Br6119"
## drop donors with limited LC; keep br 1691 (have left hemi on Xenium)
lcv <- lcv[,!(lcv$brnum %in% paste0("Br",c(5517,5276,5712)))]


### we need logcounts on the Visium object to perform NMF
lcv <- scuttle::computeLibraryFactors(lcv)
lcv <- scuttle::logNormCounts(lcv)

## load Xenium data and make sure we only have matching donors in the two datasets
lcx <- readRDS("processed-data/02_QC/02b-spe_filtered_genetargonly_log-and-nonlog-norms.RDS")
unique(lcx$brnum) # donor WITH hemisphere
lcx$donor <- lcx$brnum
lcx$donor[lcx$donor=="Br2305L"] <- "Br2305"
lcx$donor[lcx$donor=="Br1039L"] <- "Br1039"
lcx$donor[lcx$donor=="Br5854L"] <- "Br5854"
lcx$donor[lcx$donor=="Br1691L"] <- "Br1691"

stopifnot(all(unique(lcx$donor) %in% unique(lcv$brnum)))
stopifnot(all(unique(lcv$brnum) %in% unique(lcx$donor)))

## isoforms of XBP1 are separate probes for xenium, so we want to add a row that just aggregates the two. 
## we recalculate logcounts during the per-sample loop below anyhow 
## so we can get rid of the normalized counts in the xenium SPE too.
assays(lcx)[[3]] <- NULL
assays(lcx)[[2]] <- NULL
gc(full=T)


## add a row for XBP1 gene level to rowData and counts, remove isoforms
newg <- as.data.frame(rowData(lcx))
newg[(nrow(newg)+1),"Symbol"] <- "XBP1"
rownames(newg)[nrow(newg)] <- "XBP1"
newg <- as.data.table(newg,keep.rownames=T)
newg[rn=="XBP1",ID:="ENSG00000100219"]
newg[rn=="XBP1",Type:="Gene Expression"]
set(newg,i=nrow(newg),j = c(5:12),value=FALSE)
newg <- newg[!(rn %in% c("XBP1u","XBP1s"))]
newg <- DataFrame(newg,row.names=newg$rn)
newg$rn <- NULL

## remove isoforms from counts matrix
newr <- t(as.matrix(colSums(counts(lcx)[c("XBP1s","XBP1u"),],na.rm=T)))
newct <- as(rbind(as.matrix(counts(lcx)),newr),"dgCMatrix")
rownames(newct)[nrow(newct)] <- "XBP1"
newct <- newct[!(rownames(newct) %in% c("XBP1u","XBP1s")),]

## make sure we have the same genes in both of these matrices, in the same order by row and col
stopifnot(all(rownames(newct)==rownames(newg)))
stopifnot(all(colnames(newct)==colnames(lcx)))

lcx2 <- SpatialExperiment(colData=colData(lcx),rowData=newg,assays=list(counts=newct),spatialCoords = spatialCoords(lcx))

## we need the Symbol column to be titled gene_name for spaTransfer
colnames(rowData(lcx2))[2] <- "gene_name"

rm(lcx,newct,newg,newr)
gc(full=T)

## nmf projection is done using the entire visium dataset per xenium SAMPLE (encoded in brnum, 
## where one tissue section is a unique donor-hemisphere in this dataset). 
## using Cindy's implementation of this from https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/devel/code/analysis/04_label_transfer/01_run_spaTransfer.R

## the alpha parameter (elasticnet mixing parameter?) needs to be specified for this to work.
## based on the example scripts in the draft repo for the spatransfer manuscript and per cindy, use alpha=0. 

## split xenium into per-sample SPE
lcxes <- lapply(unique(lcx2$brnum),function(x){
    lcx2[,lcx2$brnum==x]
})
names(lcxes) <- unique(lcx2$brnum)
lcxes <- lapply(lcxes,scuttle::computeLibraryFactors)
lcxes <- lapply(lcxes,scuttle::logNormCounts)

rm(lcx2)
gc(full=T)

## technicalVarName should be the sample identifier in the source(reference) dataset.
## use tol 1e-8 (vs default 1e-5) for higher stringency (passed to nmfLabelTransfer::run_nmf)
## use k=63, the crossVal-checked optimal number of nmf factors

# use 15 threads; prevent data.table parallelization overcommiting threads using setDTthreads
setDTthreads(1,restore_after_fork=FALSE) # prevent resource competition
options(RcppML.verbose = TRUE)
options(RcppML.threads=15)
vis_to_xen <- transfer_labels(targets = lcxes,
	source = lcv,
	assay = "logcounts",
	annotationsName = "label",
	technicalVarName = "sample_id",
	seed = 42,
	tol=1e-8,
	k=63,
	threads=15,
	alpha=0)

saveRDS(vis_to_xen,"processed-data/99_UsingVisSensitivityDomains/01_labelTransfer/05b-LCdonsonly_VistoXen_NMFlabelxfer_tol1e8_k63_alph0.RDS")

## sessionInf
sessionInfo()

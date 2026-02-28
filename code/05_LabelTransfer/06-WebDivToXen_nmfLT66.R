setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC')

# Load libraries
library(data.table)
library(Matrix) # otherwise rcppml gives error 'no item called "package:Matrix" on the search list'
library(RcppML)
library(SingleCellExperiment)
library(SpatialExperiment)
library(nmfLabelTransfer)

lcsn <- readRDS("processed-data/05_LabelTransfer/WebDivLCSN_copy_for_labelxfer.RDS")

### we already have logcounts here

## load Xenium data and make sure we only have matching donors in the two datasets
lcx <- readRDS("processed-data/02_QC/02b-spe_filtered_genetargonly_log-and-nonlog-norms.RDS")

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

## nmf projection is done using the entire lcsn dataset per xenium SAMPLE (encoded in brnum, 
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
## use k=66, the crossVal-checked optimal number of nmf factors

# use 15 threads; prevent data.table parallelization overcommiting threads using setDTthreads
setDTthreads(1,restore_after_fork=FALSE) # prevent resource competition
options(RcppML.verbose = TRUE)
options(RcppML.threads=15)
lcsn_to_xen <- transfer_labels(targets = lcxes,
	source = lcsn,
	assay = "logcounts",
	annotationsName = "label",
	technicalVarName = "Sample",
	seed = 42,
	tol=1e-8,
	k=66,
	threads=15,
	alpha=0)

saveRDS(lcsn_to_xen,"processed-data/05_LabelTransfer/06-WebDivToXen_NMFlabelxfer_tol1e8_k66_alph0.RDS")

## sessionInf
sessionInfo()

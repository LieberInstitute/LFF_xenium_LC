#
#
#
#
#
#
#
knitr::opts_chunk$set(fig.height = 10,fig.width = 7,include = FALSE)
#### sets tab autocompletion for directories to begin from the directory containing the .Rproj session, instead of the script's directory if different
knitr::opts_knit$set(root.dir = here::here())

library(data.table)
library(ggplot2)
library(ggtext)
library(gridExtra)
library(SpatialExperiment)
library(ggrastr)
library(escheR)
library(jaffelab)
# library(SpatialFeatureExperiment)

## rstudio GUI tweaks
options('styler.addins_style_transformer' = 'biocthis::bioc_style()')
##

## >= and <= comparators that account for floating point error
greaterOrEqual <- function(x,y) {
precision <- sqrt(.Machine$double.eps)
(x >= y) | (abs(x-y) <= precision)
}

lessOrEqual <- function(x,y){
precision <- sqrt(.Machine$double.eps)
(x <= y) | (abs(x-y) <= precision)}

## enable forked parallel processing with BiocParallel::multicoreParam, future::, etc. seemed to need this a couple times, but otherwise havent so its here as a preventative measure. part of this is adding the line 
# OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES
# to Renviron.site. see e.g. top response on https://stackoverflow.com/questions/73638290/python-on-mac-is-it-safe-to-set-objc-disable-initialize-fork-safety-yes-globall 
library(parallelly)
options(parallelly.supportsMulticore.disableOn='')
options(parallelly.fork.enable=TRUE)
library(BiocParallel)
options(bphost='localhost')


## ggplot defaults
theme_set(theme_bw()+theme(axis.text.x = element_text(size = 8,color = '#000000'), axis.title.x = element_text(size = 9,color = '#000000'), axis.text.y = element_text(size = 8,color = '#000000'), axis.title.y = element_text(size =9,color = '#000000'), plot.title = element_markdown(size = 10,hjust=0.5,color = '#000000'), strip.text = element_text(size=11), legend.text = element_text(size=7,color = '#000000'), legend.title = element_text(size=8,hjust=0.5,color = '#000000'),axis.ticks = element_line(linewidth=0.5),panel.grid.major = element_line(linewidth=0.5),panel.grid.minor=element_line(linewidth=0.25),plot.title.position='plot'))

## last of all, unload the base package datasets, whose data keep getting in the way of autocompletions (e.g., Theoph is priortized over TRUE)
unloadNamespace('datasets')
#
#
#
#
#
#
#
pal <- readRDS("../LFF_spatial_LC/plots/finalpal_somuchfornorainbows.RDS")

xlc <- readRDS("processed-data/02_QC/02b-spe_filtered_genetargonly_log-and-nonlog-norms.RDS")

lxf <- fread("processed-data/05_LabelTransfer/05b-LCsmoothings_LUT.txt")
## use the strict domain for plotting
lxf <- lxf[,.(rn,domain12252,nmfpred)]
## convert working names to plotting names for cluster labels
lxf[nmfpred=="Astro_Oligo",nmfpred:="Astro-Oligo Mixed"]
lxf[nmfpred=="Oligo_Astro",nmfpred:="Oligo-Astro Mixed"]
lxf[nmfpred=="ACHE_SERT_Npep",nmfpred:="*ACHE*, *SERT*,<br>Neuropept. Mixed"]
lxf[nmfpred=="ACHE_SERT_LC_WM",nmfpred:="*ACHE*, *SERT*,<br>LC, Oligo Mixed"]
lxf[nmfpred=="Oligo_Complex",nmfpred:="Oligo-Other<br>Types Mixed"]

names(pal)[2] <- "*ACHE*, *SERT*,<br>LC, Oligo Mixed"
names(pal)[3] <- "*ACHE*, *SERT*,<br>Neuropept. Mixed"
names(pal)[9] <- "Oligo-Other<br>Types Mixed"

stopifnot(all(lxf$nmfpred %in% names(pal)))

## add to coldata
tmpcd <- as.data.table(colData(xlc),keep.rownames=T)
tmpcd <- merge.data.table(tmpcd,lxf,by="rn")
tmpcd <- DataFrame(tmpcd,row.names=tmpcd$rn)[colnames(xlc),]
tmpcd$rn <- NULL
colData(xlc) <- tmpcd

rm(lxf)
#
#
#
#
nmfpredpans <- list(Br5941=xlc[,xlc$brnum=="Br5941"],Br6538=xlc[,xlc$brnum=="Br6538"])

sizepal <- c(0.5,rep(0.2,8))
names(sizepal) <- names(pal)

lapply(names(nmfpredpans),function(n){
    x <- nmfpredpans[[n]]
    plt <- make_escheR(x,y_reverse=FALSE)
    plt <- plt |> add_fill("nmfpred",size=0.2,point_size = 0.2)
   
    ## override the escheR size settings to pass our size pal
    plt@layers[[1]]$aes_params$size <- NULL
    plt <- plt+aes(size=nmfpred)+
      scale_size_manual(values=sizepal)
    
    p <- plt+   #guides(color="none")+
    labs(fill="Transferred\nVisium Label")+
    guides(fill=guide_legend(override.aes=list(size=1.5)),size="none")+
    scale_fill_manual(values=pal)
    if (n=="Br6538"){
        p <- p + ggtitle("Br6538 (E2/E3), EA, M)")
    } else if (n=="Br3974"){
        p <- p + ggtitle("Br3974 (E3/E4), AA, M)")
    }
    p <- p+
    theme(panel.grid = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_markdown(size=6,margin=margin(0,0,0,0,"in")),
        legend.key.size = ggplot2::unit(0.075,"in"),
        legend.title=element_text(size=7,hjust=0.5),
        legend.background = element_blank(),
        legend.box = element_blank(),
        legend.box.spacing = ggplot2::unit(0,"in"),
        legend.key = element_blank(),
        legend.key.height=ggplot2::unit(0.0125,"in"),
        legend.key.width=ggplot2::unit(0.0125,"in"),
        legend.key.spacing.y = ggplot2::unit(0.04,"in"),
        title=element_text(size=10,hjust=0.5))
        # legend.position.inside = c(0.19,0.165),
        # legend.position="inside")
   p <- rasterize(p,layers="Points",dpi=600,dev="cairo_png")
   pdf(paste0("plots/manuscript/Fig5_astrogenes_in_LC/Fig5W-Predlabels_",n,".pdf"),height=2,width=3.25)
   print(p)
   dev.off()
})

rm(nmfpredpans)
#
#
#
#
xactxfer <- readRDS("processed-data/06_DomainDE/07-visLabsAsXfered_voomLmFit_objs_and_DE_tables_LabSex_LabE2vE4_LabSingleGenos_LabAnces_LabAncesE4E2_withinLCVisLabs.RDS")

xact.12252 <- xactxfer$E2E4label$voomLFt$EList$E[c("DBH","PLPP3","FGFR3","NTSR2"),]
# reconfigure data for jaffelab::cleaningy
# get sample metadata in the same arrangement for LC and astro pbulks
keepsamps <- grep(colnames(xact.12252),pattern="^LC_Br|^Astro_Br",value=T)
xcd <- xactxfer$E2E4label$dgel$samples[keepsamps,]
xact.12252 <- xact.12252[,keepsamps]

# plot: domain 8266
# ensure factors are set for model matrix
xcd$APOEsimple <- factor(xcd$APOEsimple,levels = c("E2","E4"))
xcd$Sex <- factor(xcd$Sex)
xcd$Age <- as.numeric(xcd$Age)
xcd$YRI <- as.numeric(xcd$YRI)
# model matrix

modmat <- model.matrix(~APOEsimple+Sex+Age+YRI,data = xcd)
## correcting for all covars (i.e., retaining no variables in the model matrix)
testres <- jaffelab::cleaningY(xact.12252,modmat,1)


pltg <- as.data.table(cbind(t(testres),xcd$label,xcd$brnum))
setnames(pltg,c((ncol(pltg)-1):ncol(pltg)),c("label","brnum"))
pltg <- as.data.table(cbind(pltg[,lapply(.SD,as.numeric),.SDcols=c(1:(ncol(pltg)-2))],pltg[,.(label,brnum)]))

pltg <- melt.data.table(pltg,id.vars=c("label","brnum"))
pltg[,titleg:=paste0("*",variable,"*")]
pltg[,label:=as.character(label)]

pdf("plots/manuscript/Fig5_astrogenes_in_LC/Fig5X-bulk_Xpr_DBH_and_astroDEGs.pdf",height=3,width=2.5)
ggplot(pltg,aes(x=label,y=value,col=label))+
    geom_boxplot(outliers = F,linewidth=0.375)+
    geom_point(alpha=0.7,size=1,stroke = 0.1,position = position_jitter(width=0.2))+
    facet_wrap(.~titleg,nrow=2)+
    xlab("Cell Label Transferred from Visium")+
    ylab("Pseudobulk Expression (Corrected for Age,\nSex, APOE Status, Genomic Ancestry)")+
    scale_color_manual(values=pal)+
    guides(color="none")+
    theme(panel.grid.major = element_line(linewidth=0.25),
         panel.grid.minor=element_line(linewidth=0.125),
         axis.text.x=element_text(size=8),
         axis.text.y=element_text(size=8),
         axis.title.x=element_text(size=9),
         axis.title.y=element_text(size=9),
         strip.text = element_markdown(size=9,margin=margin(0,0,-0.025,0,"in")),
         axis.ticks = element_line(linewidth=0.25),
         strip.background = element_blank())
dev.off()

rm(pltg,modmat,testres,xact.12252,xcd,xactxfer)
#
#
#
#
smoove <- xlc[,xlc$brnum=="Br5941"]
dompal <- c(pal[1],"#C3C3C3")
names(dompal) <- c("LC","Other")
sizepal <- c(0.5,0.25)
names(sizepal) <- c("LC","Other")

plt <- make_escheR(smoove,y_reverse=FALSE)
plt <- plt |> add_fill("domain12252",size=0.2,point_size = 0.2)

## override the escheR size settings to pass our size pal
plt@layers[[1]]$aes_params$size <- NULL
plt <- plt+aes(size=domain12252)+
    scale_size_manual(values=sizepal)

p <- plt+ 
labs(fill="Domain")+
ggtitle("Br5941 (E4/E4\nAA, M)")+
guides(fill=guide_legend(override.aes=list(size=1.5)),size="none")+
scale_fill_manual(values=dompal)+
theme(panel.grid = element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks = element_blank(),
    legend.text = element_markdown(size=6,margin=margin(0,0,0,0,"in")),
    legend.key.size = ggplot2::unit(0.075,"in"),
    legend.title=element_text(size=7,hjust=0.5),
    legend.background = element_blank(),
    legend.margin = margin(0,-0.3,0,0,"in"),
    legend.box = element_blank(),
    legend.box.spacing = ggplot2::unit(0,"in"),
    legend.key = element_blank(),
    legend.key.height=ggplot2::unit(0.0125,"in"),
    legend.key.width=ggplot2::unit(0.0125,"in"),
    legend.key.spacing.y = ggplot2::unit(0.04,"in"),
    title=element_text(size=10,hjust=0.5))
    # legend.position.inside = c(0.19,0.165),
    # legend.position="inside")
p <- rasterize(p,layers="Points",dpi=600,dev="cairo_png")
pdf(paste0("plots/manuscript/Fig5_astrogenes_in_LC/Fig5Y-SmoothedLC_Br5941.pdf"),height=2.25,width=1.5)
print(p)
dev.off()
#
#
#
#
domde <- readRDS("processed-data/06_DomainDE/02-voomLmFit_objs_and_DE_tables_LabSex_LabE2vE4_LabSingleGenos_LabAnces_LabAncesE4E2_XferedVisLabels.RDS")

dom.12252 <- domde$dom12252$E2E4label$voomLFt$EList$E[c("PLPP3","FGFR3","NTSR2"),]
# reconfigure data for jaffelab::cleaningy
# get sample metadata in the same arrangement for LC and astro pbulks
keepsamps <- grep(colnames(dom.12252),pattern="LC_Br",value=T)
xcd <- domde$dom12252$E2E4label$dgel$samples[keepsamps,]
dom.12252 <- dom.12252[,keepsamps]

# plot: domain 8266
# ensure factors are set for model matrix
xcd$APOEsimple <- factor(xcd$APOEsimple,levels = c("E2","E4"))
xcd$Sex <- factor(xcd$Sex)
xcd$Age <- as.numeric(xcd$Age)
xcd$YRI <- as.numeric(xcd$YRI)

modmat <-  model.matrix(~0+APOEsimple+Age+YRI+Sex,data = xcd)
degs <- jaffelab::cleaningY(dom.12252,modmat,2) # correct for terms that arent the two levels of APOEsimple

pltg <- as.data.table(cbind(t(degs),xcd$APOEsimple,xcd$brnum))
setnames(pltg,c((ncol(pltg)-1):ncol(pltg)),c("APOE","brnum"))
pltg <- as.data.table(cbind(pltg[,lapply(.SD,as.numeric),.SDcols=c(1:(ncol(pltg)-2))],pltg[,.(APOE,brnum)]))

pltg <- melt.data.table(pltg,id.vars=c("APOE","brnum"))
pltg[,titleg:=paste0("*",variable,"*")]
pltg[,APOE:=as.character(APOE)]
pltg[,APOE:=ifelse(APOE=="1","E2","E4")]

pdf("plots/manuscript/Fig5_astrogenes_in_LC/Fig5Z-Pbulk_astrogenes_byAPOE_inLCdomain12252.pdf",height=1.5,width=3)
ggplot(pltg,aes(x=APOE,y=value,col=APOE))+
    geom_boxplot(outliers = F,linewidth=0.375)+
    geom_point(alpha=0.7,size=1,stroke = 0.1,position = position_jitter(width=0.2))+
    facet_wrap(.~titleg,nrow=1,scales = "free_y")+
    xlab("APOE Status")+
    ylab("LC Domain Expression")+
    scale_color_manual(values=c("E2"="#009E73","E4"="#D55E00"))+
    guides(color="none")+
    theme(panel.grid.major = element_line(linewidth=0.25),
         panel.grid.minor=element_line(linewidth=0.125),
         axis.text.x=element_text(size=8),
         axis.text.y=element_text(size=8),
         axis.title.x=element_text(size=9),
         axis.title.y=element_text(size=9),
         strip.text = element_markdown(size=9,margin=margin(0,0,-0.025,0,"in")),
         axis.ticks = element_line(linewidth=0.25),
         strip.background = element_blank())
dev.off()
#
#
#
#
sessionInfo()
sessioninfo::session_info()
```
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

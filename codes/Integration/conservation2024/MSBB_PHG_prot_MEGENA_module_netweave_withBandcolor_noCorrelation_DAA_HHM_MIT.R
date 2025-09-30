##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/adgwas/hyperGeTest"
setwd(wd)
dir.create("results2024")

# set backends
colfuncBlue=colorRampPalette(c("white", "blue"))
colfuncBlue00=colorRampPalette(c("red", "blue"))
colfuncBrown=colorRampPalette(c("white", "brown"))
colfuncHeat=function(n) rev(heat.colors(n))
nCol=50 #number of colors in a gradient
#colfuncHeat(50)
coll = colfuncBlue00(50)
#
library(tidyverse)
library(reshape2)
library("NetWeaver")
library(RColorBrewer)
options(stringsAsFactors=FALSE)

setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA/FET.mod_olap")
modules <- read.table("MSBB_proteomics_MEGENA_module_concordant_in_other_omics_pres6_olap40_in_both_nets.txt",
                     sep="\t", stringsAsFactors = F, header = T)
modules$isConcordant2rosmap.prot.score = ifelse(modules$isConcordant2rosmap.prot == "Yes", -0.5,1)
modules$isConcordant2msbb.gene.score = ifelse(modules$isConcordant2msbb.gene == "Yes",-0.5,1)
modules$isConcordant2rosmap.gene.score = ifelse(modules$isConcordant2rosmap.gene == "Yes",-0.5,1)

modules <- modules[order(modules$moduleRankingOrder, decreasing = F),]
modules <- modules[1:50,]
modules$moduleRankingOrder = sqrt(1/modules$moduleRankingOrder)

#wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/adgwas/hyperGeTest"
setwd(wd)
dam = read.table("EnrichmentTable_PHG_prot_MEGENA_mod51_for_DAM_DAA.txt",
                 sep="\t", stringsAsFactors = F, header = T)
head(dam)
dam00 = dam[,c(1,2,9)]
head(dam00)
#melt(data-frame, na.rm = FALSE, value.name = “name”, id = 'columns')
df = spread(dam00, key = "Set2", value = "p.value.FDR")
head(df)

mg = read.table("EnrichmentTable_PHG_prot_MEGENA_mod50_for_MG_state.txt",
                 sep="\t", stringsAsFactors = F, header = T)
head(mg)
mg00 = mg[,c(1,2,9)]
head(mg00)
dff = spread(mg00, key = "Set2", value = "p.value.FDR")
head(dff)

df0 = merge(df, dff, by.x = 1, by.y = 1, all = F)

twas = read.table("EnrichmentTable_PHG_prot_MEGENA_mod50_for_MIT_TWAS.txt",
                sep="\t", stringsAsFactors = F, header = T)
head(twas)
twas00 = twas[,c(1,2,9)]
head(twas00)
dfff = spread(twas00, key = "Set2", value = "p.value.FDR")
head(dfff)
df1 = merge(df0, dfff, by.x = 1, by.y = 1, all = F)
head(df1)
df2 = df1[, colnames(df1)%in%c("Set1","DAA","HHM","MG2","MG3","MG5","MG10","MG11","MG12","TWAS")]
head(df2)
modules = merge(modules, df2, by.x = 1, by.y = 1, all.x = T )
head(modules)
modules = modules[order(modules$moduleRankingOrder, decreasing = T),]

# plotting
cyto=data.frame(Chr=modules$Module.id, Start=1, End=50, BandColor="black", stringsAsFactors=FALSE)
nb.cols <- 50
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
cyto$BandColor <- mycolors[1:50]

out.file <- "results2024/MSBB_protein_specific_modules_01142024_olap40_reverse_MIT_in_both_nets.pdf";
pdf(file=out.file, height=9, width=9, compress=TRUE);
rc.initialize(cyto, num.tracks=20, params=list(chr.padding=0.2, slice.size = 360)) # default 360
(params=rc.get.params())

rc.plot.area(size=0.90)

rc.plot.histogram(cyto, track.id=2, color.col='BandColor', track.border=NA, polygon.border=NA)
chrom.alias=cyto$Chr
names(chrom.alias)=cyto$Chr
rc.plot.ideogram(track.ids=1:2, plot.band=T, plot.chromosome.id=TRUE, cex.text=0.75, chrom.alias=chrom.alias,
                 track.border=NA, polygon.border=NA,color.chromosome.id = coll, las = 2)

Rank=data.frame(cyto[,c("Chr","Start","End")],
                Score= as.numeric(modules$moduleRankingOrder),
                stringsAsFactors=FALSE)
params$color.hist = "blue"
rc.plot.histogram(Rank, track.id=5, 
                  data.col="Score",
                  color.gradient = "blue",
                  fixed.height=FALSE,
                  custom.track.height=params$track.height*3,
                  track.border=NA,
                  polygon.border=NA)

# plot correlation
eigenCorCols = colnames(modules)[grep("isConcordant2",colnames(modules))] #select columns with pattern Rho
eigenCorCols = eigenCorCols[grepl("score",eigenCorCols)]
maxCorr=1 #set maximum correlation coefficient
track.border="#999999"
  track.color="white"
    track.num=6
    for(eigenCorCol in eigenCorCols){
      track.num <- track.num+1
      data.col <- 4
      color.col <- 5
      HistData=data.frame(Chr=modules$Module.id,Start=1,End=50,Data=modules[,eigenCorCol],
                          stringsAsFactors=FALSE)
      #plot positive correlation
      HistData1=HistData[HistData$Data > 0,]
      HistData1$col= "white" #colfuncBrown(nCol)[pmax(1,floor(HistData1[,data.col]*nCol/maxCorr))]
      rc.plot.histogram(HistData1, track.num, data.col, color.col=color.col,
                        fixed.height=TRUE, track.color=track.color, track.border=track.border, polygon.border=NA)
      #plot negative correlation
      HistData2=HistData[HistData$Data <= 0,]
      HistData2$Data=abs(HistData2$Data)
      HistData2$col= "#6b6b6b" #colfuncBlue(nCol)[pmax(1,floor(HistData2[,data.col]*nCol/maxCorr))]
      rc.plot.histogram(HistData2, track.num, data.col, color.col=color.col, fixed.height=TRUE,
                        track.border=track.border, polygon.border=NA)
    }
# plot heatmap
track.border="#999999"
track.color="white"
track.num=9
y.cor=rc.get.coordinates(1,1,1)$y[1]-0.05
x.cor=params$radius*0.15
bht=0.25
bwt=0.25

track.num=track.num+1
heatmapData=t(modules[,(ncol(modules)-8):ncol(modules)])
colnames(heatmapData)=modules$Module.id
#convert P value to -log10 scale
heatmapData[,]=as.integer(-log(heatmapData+1.0e-320,base=10))
#cap the maximum -log10 P value so that there is a richer color pattern
maxLogPval = 2
heatmapData[heatmapData>maxLogPval]=maxLogPval
#

rc.plot.heatmap(heatmapData, track.num, color.gradient= colfuncHeat(nCol),
                track.color=track.color, track.border=track.border, polygon.border=NA)
track.num=track.num+nrow(heatmapData)
maxLogPval = 2
rc.plot.grColLegend(x.cor+0.0, y.cor-bht, colfuncHeat(nCol), at=c(1,floor(nCol/2),nCol),
                    legend=signif(c(0,maxLogPval/2,maxLogPval),2), title=expression(paste(-log[10],"(FDR)")),
                  width=bwt, height=2*bht, cex.text=0.75)

rc.plot.track.id(3, labels=1, col="black", custom.track.height=params$track.height*2, cex=0.8)
#rc.plot.track.id(seq(7, track.num-1,by=2), labels=seq(7,track.num-1,by=2)-2, col="black", cex=0.8)
dev.off()
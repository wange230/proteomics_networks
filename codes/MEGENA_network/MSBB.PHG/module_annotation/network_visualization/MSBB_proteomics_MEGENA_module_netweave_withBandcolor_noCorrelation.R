##
rm(list=ls())
colfuncBlue=colorRampPalette(c("white", "blue"))
colfuncBlue00=colorRampPalette(c("red", "blue"))
colfuncBrown=colorRampPalette(c("white", "brown"))
colfuncHeat=function(n) rev(heat.colors(n))
nCol=50 #number of colors in a gradient
#colfuncHeat(50)
coll = colfuncBlue00(30)

library("NetWeaver")
library(RColorBrewer)
options(stringsAsFactors=FALSE)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/hyperGeTest_new/PHG_protein_mol_ranking")
modules <- read.table("PHG_proteomics_DEP_FC11_based_MEGENA_module_ranking_06_01_2020.txt",
                      sep="\t", stringsAsFactors = F, header = T)
#modules <- modules[order(modules$moduleRankScore, decreasing = T),]
modules$moduleRankWithinRegion = sqrt(1/modules$moduleRankingOrder)
modules <- modules[1:30,]
cyto=data.frame(Chr=modules$Module.id, Start=1, End=100, BandColor="black", stringsAsFactors=FALSE)
cyto=data.frame(Chr=modules$Module.id, Start=1, End=100, BandColor="blue", stringsAsFactors=FALSE)
nb.cols <- 30
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
cyto$BandColor <- mycolors[1:30]
#out.file <- "mRNA_specific_signatures.pdf";
#pdf(file=out.file, height=8, width=8, compress=TRUE);
out.file <- "MSBB_proteomics_module_ranking_withDEPs.pdf";
pdf(file=out.file, height= 5, width= 5);
rc.initialize(cyto, num.tracks=36, params=list(chr.padding=0.2, slice.size = 90)) # default 360
(params=rc.get.params())

rc.plot.area(size=1)
rc.plot.histogram(cyto, track.id=2, color.col='BandColor', track.border=NA, polygon.border=NA)
#chrom.alias=1:nrow(cyto)
#names(chrom.alias)=cyto$Chr
chrom.alias=cyto$Chr
names(chrom.alias)=cyto$Chr
rc.plot.ideogram(track.ids=1:2, plot.band=T, plot.chromosome.id=TRUE, cex.text=1.0, chrom.alias=chrom.alias,
                 track.border=NA, polygon.border=NA,color.chromosome.id = coll, las = 2)

Rank=data.frame(cyto[,c("Chr","Start","End")], Score= as.numeric(modules$moduleRankWithinRegion), stringsAsFactors=FALSE)
rc.plot.histogram(Rank,
                  track.id=5, 
                  data.col="Score", 
                  color.gradient = "blue",
                  fixed.height=FALSE,
                  custom.track.height=params$track.height*3,
                  track.border=NA,
                  polygon.border=NA )

track.border="#999999"
track.color="white"
track.num=6

y.cor=rc.get.coordinates(1,1,1)$y[1]-0.05
x.cor=params$radius*0.35
bht=0.6
bwt=0.2

track.num=track.num+1
heatmapData=t(modules[,c(6:24)])
colnames(heatmapData)=modules$Module.id
#convert P value to -log10 scale
heatmapData[,]=as.integer(-log(heatmapData+1.0e-320,base=10))
#cap the maximum -log10 P value so that there is a richer color pattern
maxLogPval=25
heatmapData[heatmapData>maxLogPval]=maxLogPval
#

rc.plot.heatmap(heatmapData, track.num, color.gradient=colfuncHeat(nCol),
                track.color=track.color, track.border=track.border, polygon.border=NA)
track.num=track.num+nrow(heatmapData)

rc.plot.grColLegend(x.cor+0.8, y.cor-bht, colfuncHeat(nCol), at=c(1,floor(nCol/2),nCol),
                    legend=signif(c(0,maxLogPval/2,maxLogPval),2), title=expression(paste(-log[10],"(P)")),
                    width=bwt, height=2*bht, cex.text=0.8)

rc.plot.track.id(4, labels=1, col="black", custom.track.height=params$track.height*2, cex=0.8)
rc.plot.track.id(seq(7, track.num-1,by=3), labels=seq(7,track.num-1,by=3)-3, col="black", cex=0.8)
dev.off()
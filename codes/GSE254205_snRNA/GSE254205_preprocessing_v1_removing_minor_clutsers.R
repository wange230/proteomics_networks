#
rm(list=ls())
wd = "/sc/arion/projects/adineto/apoe.GSE254205/seurats/"
setwd(wd)
list.files()

#library(SeuratData)
options(future.globals.maxSize= 6891289600)
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(patchwork) 
library(cowplot)
library(SeuratDisk)

load("GSE254205_snRNA_processed_v1_cellType.RData")

## annotating cell types
DefaultAssay(pbmc) = "SCT"
#ast
VlnPlot(pbmc, features = c("AQP4","GJA1","GFAP","ACSL1"), pt.size = 0.2, ncol = 3)
subclusters = unique(as.character(pbmc$CellType))
subclusters = subclusters[!subclusters%in%c("Olig1","Mic1", "Neu")]
pbmc.sub = subset(pbmc,idents = subclusters )
pdf(file = "plots/plot_AHNAK_v1_removing_minor_cluster.pdf")
print(VlnPlot(pbmc.sub, features = c("AHNAK"), pt.size = 0.2, ncol = 3))
print(FeaturePlot(pbmc.sub, feature="AHNAK", raster = F, pt.size = 0.5))
print(FeaturePlot(pbmc.sub, feature="AHNAK", raster = F, pt.size = 0.25, 
                  min.cutoff = 0.25, max.cutoff = 0.75))
dev.off()

#mic
VlnPlot(pbmc.sub, features = c("CSF1R","CD74","P2RY12"), pt.size = 0.2, ncol = 3)
pdf(file = "plots/featureplot_v1_removing_minor_cluster.pdf")
print(DimPlot(pbmc.sub, label = TRUE,raster=FALSE))
print(FeaturePlot(pbmc.sub, feature="MERTK", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="P2RY12", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="TREM2", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="INPP5D", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="CSF1R", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="AQP4", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="FLT1", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="MAG", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="PDGFRA", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="NRGN", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="GAD2", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="SLC17A7", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc.sub, feature="ACSL1", raster = F, pt.size = 0.2))
print(VlnPlot(pbmc.sub, features = c("ACSL1"), raster = F, pt.size = 0.2, ncol = 3))
dev.off()
#

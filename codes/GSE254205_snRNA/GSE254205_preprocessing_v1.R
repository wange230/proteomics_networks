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

#Convert("GSE254205_ad_raw.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc <- LoadH5Seurat("GSE254205_ad_raw.h5seurat",
                     meta.data = FALSE, misc = FALSE)
cts = pbmc@assays$RNA@counts
pbmc <- CreateSeuratObject(counts = cts, project = "apoe", 
                           min.cells = 0, min.features = 0)
pbmc

df = read.table("GSE254205_cell_sample_meta.txt",
                sep="\t",stringsAsFactors = F, header = T)
cel = colnames(pbmc@assays$RNA@counts)
length(intersect(df$index, cel))
row.names(df) = df$index
length(intersect(row.names(df), cel))

# adding sample meta data 
pbmc <- AddMetaData(object = pbmc, metadata = df)
dat = pbmc@meta.data
all.equal(dat$index, colnames(pbmc))
pbmc
# preprocessing
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-|^Mt-|^mt-") 
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)
pbmc
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt",verbose = TRUE) # %>%RunPCA()
pbmc
save(pbmc, file="GSE254205_snRNA_processed_v1.RData")

# batch correction
#library(harmony)
load("GSE254205_snRNA_processed_v1.RData")
pbmc = pbmc %>%RunPCA()
Idents(pbmc) = "sample"
#pbmc <- RunHarmony(object = pbmc,
 #                        group.by.vars = 'sample',
  #                       reduction = 'pca',
   #                      assay.use = 'SCT',
    #                     project.dim = FALSE,
     #                    reduction.save = "harmony.sct"
      #                   ,plot_convergence = TRUE)     
#pbmc <- RunUMAP(object = pbmc,
 #                     dims = 1:50,
  #                    reduction = 'harmony.sct',
   #                   reduction.name = 'umap.rna',
    #                  reduction.key = 'rnaUMAP_',
     #                 a = 1.5, b = 0.9)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc,dims = 1:30,verbose = FALSE)
pbmc <- FindClusters(pbmc, algorithm = 3, 
                     resolution = 0.075,
                     verbose = FALSE)
DimPlot(pbmc, label = TRUE,raster=FALSE)
pbmc@reductions
DimPlot(pbmc, reduction = 'umap',group.by = 'sample', 
        pt.size = 0.1,
        raster=FALSE)

# annotating cell types
# add annotations
pbmc<- RenameIdents(pbmc, '0' = 'Olig')
pbmc<- RenameIdents(pbmc, '1' = 'Ex')
pbmc<- RenameIdents(pbmc, '2' = 'Ast')
pbmc<- RenameIdents(pbmc, '3' = 'Mic')
pbmc<- RenameIdents(pbmc, '4' = 'Ex')
pbmc<- RenameIdents(pbmc, '5' = 'In')
pbmc<- RenameIdents(pbmc, '6' = 'In')
pbmc<- RenameIdents(pbmc, '7' = 'Opc')
pbmc<- RenameIdents(pbmc, '8' = 'Ex')
pbmc<- RenameIdents(pbmc, '9' = 'Endo')
pbmc<- RenameIdents(pbmc, '10' = 'Ex')
pbmc<- RenameIdents(pbmc, '11' = 'Neu')
pbmc<- RenameIdents(pbmc, '12' = 'Olig1')
pbmc<- RenameIdents(pbmc, '13' = 'Mic1')


pbmc$CellType <- Idents(pbmc)
pbmc$CellType = factor(pbmc$CellType, 
                       levels=c("Ast","Endo","Ex","In","Mic","Mic1","Neu","Olig","Olig1","Opc"))

DefaultAssay(pbmc) = "RNA"
pbmc <- NormalizeData(pbmc)
save(pbmc, file="GSE254205_snRNA_processed_v1_cellType.RData")
# load("GSE254205_snRNA_processed_v1_cellType.RData")

## annotating cell types
DefaultAssay(pbmc) = "SCT"
#ast
VlnPlot(pbmc, features = c("AQP4","GJA1","GFAP","ACSL1"), pt.size = 0.2, ncol = 3)

pdf(file = "plots/plot_AHNAK_v1.pdf")
print(VlnPlot(pbmc, features = c("AHNAK"), pt.size = 0.2, ncol = 3))
print(FeaturePlot(pbmc, feature="AHNAK", raster = F, pt.size = 0.5))
print(FeaturePlot(pbmc, feature="AHNAK", raster = F, pt.size = 0.25, 
                  min.cutoff = 0.25, max.cutoff = 0.75))
dev.off()

#mic
VlnPlot(pbmc, features = c("CSF1R","CD74","P2RY12"), pt.size = 0.2, ncol = 3)
pdf(file = "plots/featureplot_v1.pdf")
print(DimPlot(pbmc, label = TRUE,raster=FALSE))
print(FeaturePlot(pbmc, feature="MERTK", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="P2RY12", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="TREM2", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="INPP5D", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="CSF1R", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="AQP4", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="FLT1", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="MAG", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="PDGFRA", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="NRGN", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="GAD2", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="SLC17A7", raster = F, pt.size = 0.2))
print(FeaturePlot(pbmc, feature="ACSL1", raster = F, pt.size = 0.2))
print(VlnPlot(pbmc, features = c("ACSL1"), raster = F, pt.size = 0.2, ncol = 3))
dev.off()

#endo
VlnPlot(pbmc, features = c("FLT1","APOLD1", "CD34"), pt.size = 0.2, ncol = 3)


#olig
VlnPlot(pbmc, features = c("MAG","MOG","MOBP"), pt.size = 0.2, ncol = 3)

#opc
VlnPlot(pbmc, features = c("PDGFRA","VCAN","CA10"), pt.size = 0.2, ncol = 3)

#Ex
VlnPlot(pbmc, features = c("NRGN","CCK",""),pt.size = 0.2, ncol = 3) 
VlnPlot(pbmc, features = c("GRIN1","GRIN2B","SLC17A7","SLC17A6","GLUL"),
        pt.size = 0.2, ncol = 3) 
#In
VlnPlot(pbmc, features = c("GAD1","GAD2",""),pt.size = 0.2, ncol = 3)
VlnPlot(pbmc, features = c("GABBR1","GABBR2","SLC6A1","SLC32A1"),
        pt.size = 0.2, ncol = 3)

#neu
VlnPlot(pbmc, features = c("STMN2","VGF","OLMF3","SYNPR"),pt.size = 0.2, ncol = 3)
#VlnPlot(pbmc, features = c("APOE"),  pt.size = 0.2, ncol = 3) #?
VlnPlot(pbmc, features = c("DLX6-AS1",  "RELN", "CNR1",
                                 "OPRK1", "GABRB2"),  pt.size = 0.2, ncol = 3)
VlnPlot(pbmc, features = c("RAB3C", "SYT1", "KCNC2", "ATP6V1A", "ZMAT4"),
        pt.size = 0.2, ncol = 4)

VlnPlot(pbmc, features = c("RIMBP2", "CHGB", "GABRA1", "MYT1L", "PTHLH", "TAC3", "ARL4C", "TAC1"),
        pt.size = 0.2, ncol = 4)

VlnPlot(pbmc, features = c("SCG2", "ZCCHC12", "GALNTL6", "SERTM1", "VWC2L", "SNAP25", "DLX1"),
        pt.size = 0.2, ncol = 4)

VlnPlot(pbmc, features = c("GDA", "PLCXD3", "SYT4", "PGM2L1", "KIT", "SCN8A", "PCSK2", "GPR83"),
        pt.size = 0.2, ncol = 4)


## opc
VlnPlot(pbmc, features = c("HAS2", "GALR1", "NFYA", "TGFA", "MMRN1", "CRISPLD1", "XYLT1", "TNR"),
        pt.size = 0.2, ncol = 4)
VlnPlot(pbmc, features = c("RGS13", "FPGT", "TMEM64", "LRRK2", "MAMDC2"),
        pt.size = 0.2, ncol = 4)

## olig
VlnPlot(pbmc, features = c("UGT8", "PLP1", "ERMN", "CNDP1", "CLDN11", "CDH19", "TF"),
        pt.size = 0.2, ncol = 4)
VlnPlot(pbmc, features = c("FOLH1", "KLK6", "C10ORF90", "CNTN2", "SH3TC2", "ST18", "ERBB3"),
        pt.size = 0.2, ncol = 4)
VlnPlot(pbmc, features = c("MYRF", "CLCA4", "SLAIN1", "OPALIN", "CNP", "ENPP2", "CYBSR2"),
        pt.size = 0.2, ncol = 4)

#endo 
VlnPlot(pbmc, features = c("TGM2", "IFI27", "GPR116", "IFITM1", "TM4SF1"),
        pt.size = 0.2, ncol = 4)
VlnPlot(pbmc, features = c("ITIH5", "SELE", "CFH", "TM4SF18", "MECOM", "VWF", "ANXA3", "ITGA1"),
        pt.size = 0.2, ncol = 4)
#ast
VlnPlot(pbmc, features = c("GPR98", "BMPR1B", "ETNPPL", "GJB6", "FGFR3", "SLC25A18"),
        pt.size = 0.2, ncol = 4)
VlnPlot(pbmc, features = c("SLC1A2", "SDC4", "EDNRB", "ALDH1L1", "CHI3L1", "CLDN10", "AGT"),
        pt.size = 0.2, ncol = 4)

##mic
VlnPlot(pbmc, features = c("CCL3", "CCL3L3", "CCL4L1", "CCL4", "ITGAX",  "C1QB", "FOLR2"),
        pt.size = 0.2, ncol = 4)
VlnPlot(pbmc, features = c("INPP5D","TLR1", "SLA", "DHRS9", "P2RX4", "ARHGAP25", "KBTBDB", "TNFSF18", "HLA-DPA1"),
        pt.size = 0.2, ncol = 4)
#others
VlnPlot(pbmc, features = c("NEAT1"), pt.size = 0.2, ncol = 1)

## downstream analysis
# doublet finding
library(scater)
#library(loomR)
pbmc$sample_id = paste0("s",pbmc$sample)
pbmc.sce <- as.SingleCellExperiment(pbmc)
p1 <- plotExpression(pbmc.sce, features = "INPP5D", x = "ident") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- plotPCA(pbmc.sce, colour_by = "ident")
p1 + p2

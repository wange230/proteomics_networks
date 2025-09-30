##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/ahnak_proteomics"
setwd(wd)
library(stringr)
library(impute)
library(readxl)
library(limma)
library(readxl)
library(limma)
library(dplyr)
library(variancePartition)
library(doParallel)
library(foreach)
library(lme4)
library(sva)
library(Biobase)
library(DESeq2)
library(edgeR)
library(vsn)
library(ggfortify)
library(DSA)
library(reshape2)
library(ggplot2)
library(corrplot)
library(lme4)
library(limma)
library(groupCor)
library(kableExtra)
dataset <- read_excel("Penggrp_TMT_Report_AHNK_Astro_proteome_v0.1.0.xlsx",
                      skip = 3, sheet=3)
head(dataset)
df = data.frame(dataset[-1,])
head(df)

prot.info = df[,1:3]
head(prot.info)
row.names(df) = df$Protein.Accession..
datExp = df[,grepl("E4iPSC",colnames(df))]
head(datExp)

meta = data.frame(sampleID = colnames(datExp),
                  groups = gsub(".*\\.","",colnames(datExp)))
head(meta)
meta$treatment = rep(c("KD","CTRL"), each = 4)
head(meta)

#
info = meta
info$pheno = factor(info$treatment, levels=c("CTRL","KD"))
info$groups = paste0(info$pheno,".",info$groups)
mat = log2(datExp)
dim(mat)
colnames(mat) == info$sampleID # should all TRUE
matt <- as.matrix(mat) 

#normalization
colMedians <- apply(matt,2, median) 
exp01 <- sweep(matt,2,colMedians,'-')

## PCA analysis
dff = t(exp01)
tbl <- merge(info, dff, by.x = 1, by.y = 0, all = F)
row.names(tbl) = tbl$sampleID
tbl[1:5,1:6]

pca.dat = tbl[,-c(1:4)]
pca_res <- prcomp(pca.dat, scale. = TRUE)
autoplot(pca_res) 
pca_res$x
write.table(pca_res$x,file="plot/PCA_analysis.txt",sep="\t")
dir.create("plot")
#jpeg(file = "plot/AHNAK_KD_IPSC.jpeg",
 #   width = 1100, height = 1200, res = 300)
#print(autoplot(pca_res, data = tbl, colour = 'groups', shape = "pheno", size = 3))
#dev.off()

M <- cor(exp01)
setwd(wd)
#jpeg(file = "plot/AHNAK_KD_IPSC_cor_among_samples.jpeg",
 #      width = 1000, height = 1000, res = 150)
#corrplot(M, type="upper")
#dev.off()


## DEP analysis
design <- model.matrix(~0+info$pheno)
design
colnames(design) <- c("CTRL","KD")
head(design)
fit <- lmFit(exp01, design)
contrast.matrix <- makeContrasts(KD-CTRL,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# output
setwd(wd)
outputdir = "DEP"
if(dir.exists("DEP")){print("This dir has already existed!!!")}else{dir.create(outputdir)}

colnames(design)
ml<- topTable(fit2, coef = "KD - CTRL", number = Inf)
ml = merge(prot.info[,c(1,3)], ml, by.x = 1, by.y = 0, all=F )
ml = ml[order(ml$adj.P.Val),]
#write.table(ml, file= "DEP/AHNAK_KD_vs_CTRL_normalization.tsv", sep="\t", row.names = F)
# over
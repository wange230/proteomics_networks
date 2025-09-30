## GSEA of AHNAK KD DEP in the MSBB PHG proteomics MEGENA subnet of AHNAK
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/ahnak_proteomics"
setwd(wd)
dir ="GSEA"; if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
library(stringr)
library(impute)
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
library(groupCor)
library(kableExtra)
library(purrr)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
## gsea test *** testing
setwd("C:/Users/wange13/Desktop/new_literature/magma")
hg19 = read.table("NCBI37.3/NCBI37.3.gene.loc",sep="\t",stringsAsFactors = F, header = F)
head(hg19)
gen = hg19[,c(6,1)]
head(gen)

setwd(wd)
dep = read.table("DEP/AHNAK_KD_vs_CTRL_normalization.tsv",
                 sep="\t",stringsAsFactors = F, header = T)
head(dep)
dep = dep %>% distinct(GN,.keep_all = T)

setwd("C:/Users/wange13/Desktop/MSBB_new_rel/MEGENA_subnet")
net00 = read.table("PHG_BN.KDP_MEGENAsubnet_nodes_sp_Q09666_AHNK_HUMAN3layer.txt",
                 sep="\t",stringsAsFactors = F, header = T)
head(net00)
net01 = net00 %>% distinct(Symbol, .keep_all = T)
net01 = net01[,c(2,7)]
head(net01)

#
term2gene0 = merge(net01, gen, by.x =1 , by.y =1, all =F)
head(term2gene0)
term2gene00 = term2gene0[,c(2:3)]
head(term2gene00)
term2gene00 = term2gene00 %>% distinct(V1, .keep_all=T)
term2name = data.frame(term = c("Layer3"),name="AHNAK.KD")
head(term2name)


### testing
set.seed(123456)
dep01 = dep[,c(2,3,6)]
dep02 = merge(dep01,gen, by.x = 1, by.y=1,all=F )
head(dep02)
summary(abs(dep02$logFC))
dep02 = dep02[!abs(dep02$logFC)>log2(1.2) &
                abs(dep02$logFC)>summary(abs(dep02$logFC))[1],]
df <- NULL
cutoff = c(0.7)
p_cutoff = paste0("p_cutoff",cutoff)

for (nn in 1:length(cutoff)) {

dep03 = dep02[dep02$P.Value < cutoff[nn],]
dep03 = dep03[order(dep03$logFC, decreasing = T),]
head(dep03)
genes = dep03$logFC
names(genes) = dep03$V1

rel = GSEA(genes, exponent = 1, 
           nPerm = 10000, 
           minGSSize = 1, 
           maxGSSize = 5000,
           pvalueCutoff = 1,
           pAdjustMethod = "BH",
           TERM2GENE = term2gene00, 
           TERM2NAME = term2name, 
           verbose = TRUE, seed = FALSE )#,by = "DOSE")
dat = rel@result
dat$p_cutoff = rep(p_cutoff[nn], nrow(dat)) 
df = rbind(df,dat)
rm(dat)
}# nn
df$FC_cutoff.upper = rep(1.2, nrow(df))
df$FC_cutoff.lower = rep("summary(abs(dep02$logFC))[1]", nrow(df))
setwd(wd)
pdf(file="GSEA/DEP_cut07.pdf",
     width = 5.0, height = 4.5)
gseaplot2(rel,geneSetID =c(1),pvalue_table = TRUE)
dev.off()
#write.table(dep03, file="GSEA/AHNAK_DEP07_for_GSEA.txt",sep="\t", row.names = F)
#
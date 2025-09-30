#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA/FET.mod_olap"
setwd(wd)
library(readxl)
library(tidyverse)

setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/updated_10172023/v7/excel.files")
dataset <- read_excel("supplementary.data 8.xlsx",sheet=2,skip=4)
pr = data.frame(dataset)
table(pr$Zsummary.pres>2) # 0.7581522(279/368)
length(pr$module[pr$Zsummary.pres>10]) # 50
length(pr$module[pr$Zsummary.pres>10])/nrow(pr) # 0.1358696

rank = read_excel("supplementary.data 7.xlsx",sheet=4,skip=4)
rank = data.frame(rank[,c(1,2)])  
head(rank)
rank$Module.ranking.order = 1:nrow(rank)
rank30 = rank[1:30,c(1,2)]

setwd(wd)
# ROSMAP prot to MSBB prot
dat = read.table("EnrichmentTable_MEGENA_modules_ROSMAP_prot_vs_MSBB_prot.txt",
                 sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]

conserved = unique(dat$Set1[dat$p.value.FDR<0.05])
not.conserved = setdiff(unique(dat$Set1),conserved)
# conserved percent
length(conserved) # 248
length(not.conserved) # 120
length(conserved)/(length(conserved)+length(not.conserved))*100
# 67.3913

length(intersect(conserved, pr$module[pr$Zsummary.pres>2])) # 224
length(intersect(conserved, pr$module[pr$Zsummary.pres>10])) # 50


# ROSMAP prot to ROSMAP gene
setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/updated_10172023/v7/excel.files")
dataset <- read_excel("supplementary.data 8.xlsx",sheet=6,skip=4)
pr = data.frame(dataset)
table(pr$Zsummary.pres>2) # 0.5437158(199/366)
length(pr$module[pr$Zsummary.pres>10]) # 11
length(pr$module[pr$Zsummary.pres>10])/nrow(pr) # 0.03005464

setwd(wd)
dat = read.table("EnrichmentTable_MEGENA_modules_ROSMAP_prot_vs_ROSMAP_mRNA.txt",
                 sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]
head(dat)

conserved = unique(dat$Set1[dat$p.value.FDR<0.05])
not.conserved = setdiff(unique(dat$Set1),conserved)
# conserved percent
length(conserved) # 91
length(not.conserved) # 227
length(conserved)/(length(conserved)+length(not.conserved))*100
# 24.72826

length(intersect(conserved, pr$module[pr$Zsummary.pres>2])) # 74
length(intersect(conserved, pr$module[pr$Zsummary.pres>10])) # 11
#
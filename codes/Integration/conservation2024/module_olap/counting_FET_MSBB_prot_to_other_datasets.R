#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA/FET.mod_olap"
setwd(wd)
library(readxl)
library(tidyverse)

setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/updated_10172023/v7/excel.files")
dataset <- read_excel("supplementary.data 8.xlsx",sheet=1,skip=4)
pr = data.frame(dataset) 
table(pr$Zsummary.pres>2) # 0.873057(337/386)
length(pr$module[pr$Zsummary.pres>10]) # 105
length(pr$module[pr$Zsummary.pres>10])/nrow(pr) # 0.2720207

rank = read_excel("supplementary.data 4.xlsx",sheet=4,skip=4)
rank = data.frame(rank[,c(1,6,7)])  
head(rank)
#rank30 = rank[1:30,c(1,2)]
rank30 = rank[,c(1,2)]
top30 = rank30$Module.id

setwd(wd)
# MSBB prot to ROSMAP prot
dat = read.table("EnrichmentTable_MEGENA_modules_ROSMAP_prot_vs_MSBB_prot.txt",
                   sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]

conserved = unique(dat$Set2[dat$p.value.FDR<0.05])
not.conserved = setdiff(unique(dat$Set2),conserved)
# conserved percent
length(conserved) # 258
length(not.conserved) # 128
length(conserved)/(length(conserved)+length(not.conserved))*100
# 66.83938

length(intersect(conserved, pr$module[pr$Zsummary.pres>2])) # 250
length(intersect(conserved, pr$module[pr$Zsummary.pres>10])) # 103
Hypergeometric.test = conserved
module.preservation = pr$module[pr$Zsummary.pres>2]

library(ggVennDiagram)
library(ggplot2)
library("ggsci")
x <- list(`Module.preservation` = module.preservation, 
          `Hypergeometric.test` = Hypergeometric.test)
setwd(wd)
pdf(file="module_preservation_vs_Hypergeometric_test_for MSBB_prot_in_ROSMAP_prot_MEGENA.pdf")
p= ggVennDiagram(x, label_color = "red",
                 set_size = 6.5,
                 label_size = 7.5,
                 set_color = "blue",
                 label_txtWidth = 100,
                 edge_size = 1.5)
p + scale_fill_distiller(palette = "Pastel1")
dev.off()


# MSBB prot to MSBB gene
setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/updated_10172023/v7/excel.files")
dataset <- read_excel("supplementary.data 8.xlsx",sheet=4,skip=4)
pr = data.frame(dataset) 
table(pr$Zsummary.pres>2) # 0.7885117(302/383)
length(pr$module[pr$Zsummary.pres>10]) # 31
length(pr$module[pr$Zsummary.pres>10])/nrow(pr) # 0.08093995

setwd(wd)
dat = read.table("EnrichmentTable_MEGENA_modules_MSBB_prot_vs_MSBB_mRNA.txt",
                 sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]
head(dat)

conserved = unique(dat$Set1[dat$p.value.FDR<0.05])
not.conserved = setdiff(unique(dat$Set1),conserved)

# conserved percent
length(conserved) # 159
length(not.conserved) # 128
length(conserved)/(length(conserved)+length(not.conserved))*100
# 41.19171

length(intersect(conserved, pr$module[pr$Zsummary.pres>2])) # 145
length(intersect(conserved, pr$module[pr$Zsummary.pres>10])) # 29
Hypergeometric.test = conserved
module.preservation = pr$module[pr$Zsummary.pres>2]
x <- list(`Module.preservation` = module.preservation, 
          `Hypergeometric.test` = Hypergeometric.test)
#
setwd(wd)
pdf(file="module_preservation_vs_Hypergeometric_test_for MSBB_prot_in_MSBB_gene_MEGENA.pdf")
p= ggVennDiagram(x, label_color = "red",
                 set_size = 6.5,
                 label_size = 7.5,
                 set_color = "blue",
                 label_txtWidth = 100,
                 edge_size = 1.5)
p + scale_fill_distiller(palette = "Pastel1")
dev.off()
# over
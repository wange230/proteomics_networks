#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA/FET.mod_olap"
setwd(wd)
library(readxl)
library(tidyverse)

setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/updated_10172023/v7/excel.files")
dataset <- read_excel("supplementary.data 8.xlsx",sheet=1,skip=4)
pr = data.frame(dataset)
pr.sel = pr[pr$Zsummary.pres > 10,]

rank = read_excel("supplementary.data 4.xlsx",sheet=4,skip=4)
rank = data.frame(rank[,c(1,6,7)])  
head(rank)
rank30 = rank[1:30,c(1,2)]

# MSBB vs others
#top30 = c("M3","M182","M475","M750","M34","M39","M193","M492","M35","M210",
      #    "M245","M762","M46","M510","M771","M533","M9","M221","M520","M5",
 #         "M70","M769","M43","M306","M659","M806","M2","M120","M319","M186")
top30 = rank30$Module.id
tbl <- NULL
setwd(wd)
# MSBB prot to ROSMAP prot
dat = read.table("EnrichmentTable_MEGENA_modules_ROSMAP_prot_vs_MSBB_prot.txt",
                   sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]

conserved = unique(dat$Set2[dat$p.value.FDR<0.05])
not.conserved = setdiff(unique(dat$Set2),conserved)

conserved.status = ifelse(top30%in%conserved, 1, ifelse(top30%in%not.conserved,0,10))
table(conserved.status)

df = data.frame(module = top30,
                is.conserved = conserved.status,
                comparison = rep("MSBB.prot2ROSMAP.prot", length(top30))
                )
head(df)
tbl = rbind(tbl,df)
rm(df,dat)

# MSBB prot to MSBB gene
dat = read.table("EnrichmentTable_MEGENA_modules_MSBB_prot_vs_MSBB_mRNA.txt",
                 sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]
head(dat)

conserved = unique(dat$Set1[dat$p.value.FDR<0.05])
not.conserved = setdiff(unique(dat$Set1),conserved)

conserved.status = ifelse(top30%in%conserved, 1, ifelse(top30%in%not.conserved,0,10))
table(conserved.status)

df = data.frame(module = top30,
                is.conserved = conserved.status,
                comparison = rep("MSBB.prot2MSBB.gene", length(top30))
)
head(df)
tbl = rbind(tbl,df)
rm(df,dat)

# MSBB prot to ROSMAP gene
dat = read.table("EnrichmentTable_MEGENA_modules_MSBB_prot_vs_ROSMAP_mRNA.txt",
                 sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]
head(dat)

conserved = unique(dat$Set1[dat$p.value.FDR<0.05])
not.conserved = setdiff(unique(dat$Set1),conserved)

conserved.status = ifelse(top30%in%conserved, 1, ifelse(top30%in%not.conserved,0,10))
table(conserved.status)

df = data.frame(module = top30,
                is.conserved = conserved.status,
                comparison = rep("MSBB.prot2ROSMAP.gene", length(top30))
)
head(df)
tbl = rbind(tbl,df)
rm(df,dat)

head(tbl)

tbl00 = spread(tbl, key = "comparison", value = "is.conserved")
row.names(tbl00) = tbl00$module
head(tbl00)
colnames(tbl00)[2:4] =  gsub("MSBB.prot2","",colnames(tbl00)[2:4])
head(tbl00)
tbl01 = merge(rank30, tbl00, by.x = 1, by.y = 1, all = F)

row.names(tbl01) = tbl01$Module.id
mat = tbl00[,-c(1,2)]
keep = rowSums(mat) > 0
table(keep)
mat00 = mat[keep,]

# 
library(GOplot)
chord00 = tbl01[tbl01$Module.id%in%row.names(mat00),]
chord00$logFC = chord00$Module.rankingScore
chord00 = data.frame(chord00)
chord01 = chord00[,-c(1:2)]
GOChord(chord01, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)

setwd(wd)
pdf(file="summary_module_preservation_MSBB_prot_to_other_datasets_top30.pdf")
print(GOChord(chord01, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5))
dev.off()
#circ <- circle_dat(EC$david, EC$genelist)
#chord <- chord_dat(data = circ, genes = EC$genes, process = EC$process)

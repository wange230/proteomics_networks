#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA/FET.mod_olap"
setwd(wd)
library(readxl)
library(tidyverse)

setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/updated_10172023/v7/excel.files")
dataset <- read_excel("supplementary.data 8.xlsx",sheet=1,skip=4)
pr = data.frame(dataset) 
pr = pr[c("module", "module.size","Zsummary.pres")]
head(pr)

rank = read_excel("supplementary.data 4.xlsx",sheet=4,skip=4)
rank = data.frame(rank[,c(1,6,7)])  
head(rank)

msbb.prot = merge(rank, pr, by.x = 1, by.y = 1, all = F)
head(msbb.prot)
colnames(msbb.prot)[ncol(msbb.prot)] = "msbb.prot2rosmap.prot.preserved.score"
rm(rank)

setwd(wd)
# MSBB prot to ROSMAP prot
dat = read.table("EnrichmentTable_MEGENA_modules_ROSMAP_prot_vs_MSBB_prot.txt",
                   sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]
head(dat)

dat$msbb.prot.pcnt = dat$Overlap.size/dat$Sampling.size*100
dat$rosmap.prot.pcnt = dat$Overlap.size/dat$Positive.size*100
head(dat)
datt = dat[c("Set2","p.value.FDR","msbb.prot.pcnt","rosmap.prot.pcnt")]
head(datt)
dattt = datt[datt$p.value.FDR < 0.05 & (datt$msbb.prot.pcnt > 25 | datt$rosmap.prot.pcnt>25),]
splt = split(dattt, dattt$Set2)
conserved = names(splt)

msbb.prot$msbb.prot2rosmap.prot.olap.test = ifelse(msbb.prot$Module.id%in%conserved,"significant","no")
head(msbb.prot)

msbb.prot$isConcordant2rosmap.prot = ifelse(msbb.prot$msbb.prot2rosmap.prot.preserved.score > 6 &
                                              msbb.prot$msbb.prot2rosmap.prot.olap.test == "significant","Yes","No")
table(msbb.prot$isConcordant2rosmap.prot)
table(msbb.prot$isConcordant2rosmap.prot[msbb.prot$moduleRankingOrder<101])
table(msbb.prot$isConcordant2rosmap.prot[msbb.prot$moduleRankingOrder<30])
msbb.prot$Module.id[msbb.prot$isConcordant2rosmap.prot=="Yes"&msbb.prot$moduleRankingOrder<101]
msbb.prot$Module.id[msbb.prot$isConcordant2rosmap.prot=="Yes"&msbb.prot$moduleRankingOrder<31]
head(msbb.prot)
rm(pr, dat,datt,dattt,dataset,splt,conserved)

# MSBB prot to MSBB gene
setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/updated_10172023/v7/excel.files")
dataset <- read_excel("supplementary.data 8.xlsx",sheet=4,skip=4)
pr = data.frame(dataset) 
pr = pr[c("module.id", "Zsummary.pres")]
head(pr)

head(msbb.prot)
msbb.prot = merge(msbb.prot, pr, by.x = 1, by.y = 1, all.x = T)
table(is.na(msbb.prot$Zsummary.pres))
msbb.prot$Zsummary.pres = ifelse(is.na(msbb.prot$Zsummary.pres), 0, msbb.prot$Zsummary.pres)
colnames(msbb.prot)[ncol(msbb.prot)] = "msbb.prot2msbb.gene.preserved.score"

setwd(wd)
# MSBB prot to MSBB gene
dat = read.table("EnrichmentTable_MEGENA_modules_MSBB_prot_vs_MSBB_mRNA.txt",
                 sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]
head(dat)
dat$msbb.prot.pcnt = dat$Overlap.size/dat$Positive.size*100
dat$msbb.gene.pcnt = dat$Overlap.size/dat$Sampling.size*100
head(dat)
datt = dat[c("Set1","p.value.FDR","msbb.prot.pcnt","msbb.gene.pcnt")]
head(datt)
dattt = datt[datt$p.value.FDR < 0.05 & (datt$msbb.prot.pcnt > 25 | datt$msbb.gene.pcnt>25),]
splt = split(dattt, dattt$Set1)
conserved = names(splt)

msbb.prot$msbb.prot2msbb.gene.olap.test = ifelse(msbb.prot$Module.id%in%conserved,"significant","no")
head(msbb.prot)

msbb.prot$isConcordant2msbb.gene = ifelse(msbb.prot$msbb.prot2msbb.gene.preserved.score > 6 &
                                              msbb.prot$msbb.prot2msbb.gene.olap.test == "significant","Yes","No")
table(msbb.prot$isConcordant2msbb.gene)
table(msbb.prot$isConcordant2msbb.gene[msbb.prot$moduleRankingOrder<101])
table(msbb.prot$isConcordant2msbb.gene[msbb.prot$moduleRankingOrder<30])
msbb.prot$Module.id[msbb.prot$isConcordant2msbb.gene=="Yes"&msbb.prot$moduleRankingOrder<101]
msbb.prot$Module.id[msbb.prot$isConcordant2msbb.gene=="Yes"&msbb.prot$moduleRankingOrder<31]
head(msbb.prot)
rm(pr, dat,datt,dattt,dataset,splt,conserved)


# MSBB prot to ROSMAP gene
setwd("C:/Users/wange13/Desktop/MEGENA_module_conservation/MSBB_protein_to_ROSMAP_mRNA_preservation")
dataset = read.table("MEGENA_module_preservation_PHG_protein_ref_ROSMAP_mRNA_valid_perm200.txt",
                     sep="\t", stringsAsFactors = F,header = T)
pr = data.frame(dataset) 
pr = pr[c("module.id", "Zsummary.pres")]
head(pr)

head(msbb.prot)
msbb.prot = merge(msbb.prot, pr, by.x = 1, by.y = 1, all.x = T)
table(is.na(msbb.prot$Zsummary.pres))
msbb.prot$Zsummary.pres = ifelse(is.na(msbb.prot$Zsummary.pres), 0, msbb.prot$Zsummary.pres)
colnames(msbb.prot)[ncol(msbb.prot)] = "msbb.prot2rosmap.gene.preserved.score"

setwd(wd)
dat = read.table("EnrichmentTable_MEGENA_modules_MSBB_prot_vs_ROSMAP_mRNA.txt",
                 sep="\t",stringsAsFactors = F, header = T)
dat = dat[,c(1:5,9)]
head(dat)
dat$msbb.prot.pcnt = dat$Overlap.size/dat$Positive.size*100
dat$rosmap.gene.pcnt = dat$Overlap.size/dat$Sampling.size*100
head(dat)
datt = dat[c("Set1","p.value.FDR","msbb.prot.pcnt","rosmap.gene.pcnt")]
head(datt)
dattt = datt[datt$p.value.FDR < 0.05 & (datt$msbb.prot.pcnt > 25 | datt$rosmap.gene.pcnt>25),]
splt = split(dattt, dattt$Set1)
conserved = names(splt)

msbb.prot$msbb.prot2rosmap.gene.olap.test = ifelse(msbb.prot$Module.id%in%conserved,"significant","no")
head(msbb.prot)

msbb.prot$isConcordant2rosmap.gene = ifelse(msbb.prot$msbb.prot2rosmap.gene.preserved.score > 6 &
                                            msbb.prot$msbb.prot2rosmap.gene.olap.test == "significant","Yes","No")
table(msbb.prot$isConcordant2rosmap.gene)
table(msbb.prot$isConcordant2rosmap.gene[msbb.prot$moduleRankingOrder<101])
table(msbb.prot$isConcordant2rosmap.gene[msbb.prot$moduleRankingOrder<30])
msbb.prot$Module.id[msbb.prot$isConcordant2rosmap.gene=="Yes"&msbb.prot$moduleRankingOrder<101]
msbb.prot$Module.id[msbb.prot$isConcordant2rosmap.gene=="Yes"&msbb.prot$moduleRankingOrder<31]
head(msbb.prot)
rm(pr, dat,datt,dattt,dataset,splt,conserved)

write.table(msbb.prot, file="MSBB_proteomics_MEGENA_module_concordant_in_other_omics_pres6_olap25_in_both_nets.txt",
sep="\t", row.names = F, quote = F)
# over

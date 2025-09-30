##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/volcano_plot"
setwd(wd)
library(EnhancedVolcano)
library(airway)
library(magrittr)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019")
meta = read.table("protein.meta.with.uniqueSymbol.txt", sep="\t", stringsAsFactors = F, header = T)
head(meta)


cdr <- read.table("Protein_exprCorrected_HvsL_CDR_final.tsv", 
                  sep="\t", stringsAsFactors = F, header = T)
head(cdr)
cdr_meta <- merge(meta,cdr, by.x = 0, by.y =0, all =FALSE)
cdr_meta <- cdr_meta[order(cdr_meta$adj.P.Val, decreasing = F),]
row.names(cdr_meta) <- cdr_meta$protein.gene
head(cdr_meta)
dep.cdr = cdr_meta[abs(cdr_meta$logFC)>log2(1.1)&cdr_meta$adj.P.Val<0.05,]
dep.cdr$sign = ifelse(dep.cdr$logFC>0, "Up","Dn")


cerad <- read.table("Protein_exprCorrected_HvsL_CERAD_final.tsv", 
                    sep="\t", stringsAsFactors = F, header = T)
head(cerad)
cerad_meta <- merge(meta,cerad, by.x = 0, by.y =0, all =FALSE)
cerad_meta <- cerad_meta[order(cerad_meta$adj.P.Val, decreasing = F),]
row.names(cerad_meta) <- cerad_meta$protein.gene
head(cerad_meta)
dep.cerad = cerad_meta[abs(cerad_meta$logFC)>log2(1.1)&cerad_meta$adj.P.Val<0.05,]
dep.cerad$sign = ifelse(dep.cerad$logFC>0, "Up","Dn")

sig1 = dep.cdr$Protein[dep.cdr$sign=="Dn"]
sig2 = dep.cdr$Protein[dep.cdr$sign=="Up"]
sig3 = dep.cerad$Protein[dep.cerad$sign=="Dn"]
sig4 = dep.cerad$Protein[dep.cerad$sign=="Up"]
list.files("C:/pc.backup/Documents/R_code/R_functions_for_cor_network")
source("C:/pc.backup/Documents/R_code/R_functions_for_cor_network/Function_4signature_venn_fn_legend_pdf.R")
setwd(wd)

venn4sig(sig1, sig2, sig3,sig4, "CDR.Down","CDR.Up","CERAD.Down","CERAD.Up", "olap_DEPs_from_CDR_and_CERAD")




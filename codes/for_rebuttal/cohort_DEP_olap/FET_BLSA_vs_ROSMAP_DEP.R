## BLSA PFC 47 samples vs ROSMAP
rm(list=ls())
wd = "C:/Users/wange13/Desktop/DEP.validation/"
library(purrr)
library(tidyverse)

# determine background
setwd("C:/Users/wange13/Desktop/Proteomics_lita/jh_proteomics/DEP2024.final")

blsa <- read.table("PFC47samples_clean_data_normalized_LMM_adj_ADvsNL_unique_prot.tsv",
                   sep="\t", stringsAsFactors = F, header = T)
head(blsa)
blsa.gn <- unique(blsa$GN)
blsa.gn <- blsa.gn[!is.na(blsa.gn)]

blsa.dep <- blsa[blsa$adj.P.Val < 0.05, ] 

setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/output")
rosmap = read.table("ROSMAP_TNT_proteomics_expr_sex_aod_adj_lm_401.tsv",
                    sep="\t",stringsAsFactors = F, header = T)
prot = rosmap$gene.id
splt = strsplit(prot, "\\|")
ros.gn  = map_chr(splt, function(x){x[[1]]})
ros.gn = unique(ros.gn)

blsa.ros.bkg = intersect(blsa.gn,ros.gn)
length(blsa.ros.bkg) #3540

#
blsa.dep$module = ifelse(blsa.dep$logFC>0, "Up", "Down")
head(blsa.dep)
set1 <- blsa.dep[blsa.dep$GN%in%blsa.ros.bkg,c(3,ncol(blsa.dep))]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP")
all <- read.table("MSBB_ROSMAP_protein_signatures_at_protein_GN.level.txt",
                  sep="\t", stringsAsFactors = F, header = T)
ros <- all[all$cohort == "ROSMAP",]
set2 <- ros[ros$protein.id%in%blsa.ros.bkg,c(1,2)]
head(set2)
splt2 <- split(set2, set2$dn_up)
length(splt2)

source("C:/pc.backup/common_tools_datasets/IGAP_signif_p005/Sig_seperate/astrocyte/geneSetOverLap.R")
size =  length(blsa.ros.bkg) #3540

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
#write.table(enTable, file = "hyperGeTest/EnrichmentTable_BLSA_ROSMAP_union_unique_DEPs_by_GN.txt",
 #            sep="\t", row.names = FALSE)

#
## ven diagrams
ROSMAP.dep = ros
head(ROSMAP.dep)
ROSMAP.dep.neg <- unique(ROSMAP.dep[ROSMAP.dep$dn_up == "dn",]$protein.id)
ROSMAP.dep.pos <- unique(ROSMAP.dep[ROSMAP.dep$dn_up == "up",]$protein.id)

blsa.dep = blsa.dep[c("GN","module")]
blsa.dep.neg <- unique(blsa.dep[blsa.dep$module == "Down",]$GN)
blsa.dep.pos <- unique(blsa.dep[blsa.dep$module == "Up",]$GN)

setwd(wd)
dir.create("venn.plots")
wd = paste0(wd, "/","venn.plots")
source("C:/pc.backup/Documents/R_code/R_functions_for_cor_network/Function_4signature_venn_fn_legend_pdf.R")
fix(venn4sig)
venn4sig(ROSMAP.dep.neg,ROSMAP.dep.pos,blsa.dep.neg,blsa.dep.pos, 
         fn = "ROSMAPvsBLSA_DEP00","ROSMAP(Down)","ROSMAP(Up)","BLSA(Down)","BLSA(Up)")


## BLSA PFC 47 samples
rm(list=ls())
wd = "C:/Users/wange13/Desktop/DEP.validation/"

# determine background
setwd("C:/Users/wange13/Desktop/Proteomics_lita/jh_proteomics/DEP2024.final")

blsa <- read.table("PFC47samples_clean_data_normalized_LMM_adj_ADvsNL_unique_prot.tsv",
                     sep="\t", stringsAsFactors = F, header = T)
head(blsa)
blsa.gn <- unique(blsa$GN)
blsa.gn <- blsa.gn[!is.na(blsa.gn)]

blsa.dep <- blsa[blsa$adj.P.Val < 0.05, ] 

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics")
df <- read.table("MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj.tsv",
                 sep="\t", stringsAsFactors = F, header = T, quote = '"')
head(df)
msbb.gn = unique(df$Symbol)

blsa.msbb.bkg = intersect(blsa.gn,msbb.gn)
length(blsa.msbb.bkg)

#
blsa.dep$module = ifelse(blsa.dep$logFC>0, "Up", "Down")
head(blsa.dep)
set1 <- blsa.dep[blsa.dep$GN%in%blsa.msbb.bkg,c(3,ncol(blsa.dep))]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP")
all <- read.table("MSBB_ROSMAP_protein_signatures_at_protein_GN.level.txt",
                   sep="\t", stringsAsFactors = F, header = T)
phg <- all[all$cohort == "MSBB",]
set2 <- phg[phg$protein.id%in%blsa.msbb.bkg,c(1,2)]
head(set2)
splt2 <- split(set2, set2$dn_up)
length(splt2)

source("C:/pc.backup/common_tools_datasets/IGAP_signif_p005/Sig_seperate/astrocyte/geneSetOverLap.R")
size =  length(blsa.msbb.bkg) #3505

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
#write.table(enTable, file = "hyperGeTest/EnrichmentTable_BLSA_MSBB_union_unique_DEPs_by_GN.txt",
 #            sep="\t", row.names = FALSE)
#

## ven diagrams
msbb.dep = phg
head(msbb.dep)
msbb.dep.neg <- unique(msbb.dep[msbb.dep$dn_up == "dn",]$protein.id)
msbb.dep.pos <- unique(msbb.dep[msbb.dep$dn_up == "up",]$protein.id)

blsa.dep = blsa.dep[c("GN","module")]
blsa.dep.neg <- unique(blsa.dep[blsa.dep$module == "Down",]$GN)
blsa.dep.pos <- unique(blsa.dep[blsa.dep$module == "Up",]$GN)

setwd(wd)
dir.create("venn.plots")
wd = paste0(wd, "/","venn.plots")
source("C:/pc.backup/Documents/R_code/R_functions_for_cor_network/Function_4signature_venn_fn_legend_pdf.R")
fix(venn4sig)
venn4sig(msbb.dep.neg,msbb.dep.pos,blsa.dep.neg,blsa.dep.pos, 
         fn = "MSBBvsBLSA_DEP00","MSBB(Down)","MSBB(Up)","BLSA(Down)","BLSA(Up)")

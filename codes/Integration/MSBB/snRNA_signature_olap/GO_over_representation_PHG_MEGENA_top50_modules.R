#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/adgwas/hyperGeTest"
setwd(wd)
library(clusterProfiler)
library(org.Hs.eg.db)
dir.create("go.megenatop50")

df = read.table("EnrichmentTable_PHG_prot_MEGENA_mod50_for_MG_state.txt",
                 sep="\t",stringsAsFactors = F, header = T)
mods = unique(df$Set1)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein")
megena = read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                     sep="\t",stringsAsFactors = F,header = T)
head(megena)
megena = megena[,c(2,4)]
bkg = unique(megena$Symbol)
bkg = bkg[!is.na(bkg)]

setwd(wd)
for (ii in 1:length(mods)){
olap.m3 = megena$Symbol[megena$module%in%mods[ii]]
m3 = unique(olap.m3)
m3 = m3[!is.na(m3)]
if(length(m3) < 10) {next}
ego <- enrichGO(gene          = m3,
                universe      = bkg,
                keyType = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = TRUE)
tbl = ego@result
tbl = tbl[tbl$p.adjust<0.05,]
if (nrow(tbl)<1) {next}
write.table(tbl, file= paste0("go.megenatop50/","enrichrGO_",mods[ii],".txt"),
                              sep ="\t", row.names = F, quote = F)
rm(ego, tbl, m3)

} # ii

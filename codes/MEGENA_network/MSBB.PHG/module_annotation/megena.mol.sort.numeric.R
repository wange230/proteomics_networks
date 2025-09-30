##
rm(list=ls())
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein")

megena <- read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_module_summary.tsv",
                     sep="\t", stringsAsFactors = F, header = T)
head(megena)
megena_ord <- megena[order(as.numeric(gsub("M","",megena$mod.parent)), as.numeric(gsub("M","",megena$module))),]

write.table(megena_ord, file="MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_module_summary_numeric.tsv",
            sep="\t", row.names = F, quote = F)


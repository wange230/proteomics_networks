##
rm(list =ls())
wd="C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/hyperGeTest_new/PHG_protein_mol_ranking"
setwd(wd)
rank <- read.table("PHG_FC11_DEP_based_protein_MEGENA_moldule_ranking.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(rank)

big.table <- read.table("PHG_FC11_DEP_vs_MEGENA_moldule_olap_adjP.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
colnames(big.table)[2:20] <- c("HvsL_Braak.FDR.negative","HvsL_Braak.FDR.positive","HvsL_CDR.FDR.negative","HvsL_CDR.FDR.positive",
                        "HvsL_CERAD.FDR.negative","HvsL_CERAD.FDR.positive","HvsL_PlaqueMean.FDR.negative","HvsL_PlaqueMean.FDR.positive",
                        "HvsM_Braak.FDR.negative",	"HvsM_Braak.FDR.positive",	"HvsM_CDR.FDR.negative","HvsM_CDR.FDR.positive",
                        "HvsM_CERAD.FDR.negative",	"HvsM_CERAD.FDR.positive",	"HvsM_PlaqueMean.FDR.negative","HvsM_PlaqueMean.FDR.positive",
                         "MvsL_CERAD.FDR.negative",	"MvsL_CERAD.FDR.positive","MvsL_PlaqueMean.FDR.positive")

head(big.table)


setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein/MSigDB")
go <- read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_Ontology.xls", 
                 sep="\t", stringsAsFactors = FALSE, header = TRUE)
names(go)
go <- go[,c(1:4,10)]


df <- NULL
splt <- split(go, go$module)
for(nn in 1:length(splt)) {
  dat <- splt[[nn]][which(splt[[nn]]$Corrected_P == min(splt[[nn]]$Corrected_P)),]
  if(dim(dat)[1] >1) {dat <- dat[1,]}
df <- rbind(df, dat)
}
head(df)

combs <- merge(rank, df, by.x = 1, by.y = 1, all.x = TRUE)

combs00 <- merge(big.table, df, by.x = 1, by.y = 1, all.x = TRUE)

big.table <- merge(rank, combs00, by.x = 1, by.y = 1, all.x = TRUE)
big.table00 <- big.table[, c(1,25:28, 6:24, 2:5)] 

setwd(wd)
write.table(big.table00, file = "PHG_proteomics_DEP_FC11_based_MEGENA_module_ranking_02_20_2020_updated.txt",
            sep="\t", row.names = FALSE, quote = F)





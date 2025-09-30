##
rm(list=ls())
wd= "C:/Users/wange13/Desktop/MSBB_new_rel/hyperGeTest2022"
setwd(wd)

# adding HvsM
mol <- read.table("EnrichmentTable_PHG_proteomics_bnGlobalKDA_nodes_sig_over_DEP_HvsL.txt", 
                      sep="\t", stringsAsFactors = FALSE, header = TRUE)
mol<- mol[,c(1,2,9)]
head(mol)
splt <- split(mol, mol$Set1)
df = data.frame(KDP = unique(mol$Set2))
for (nn in 1:length(splt)) {
dat = splt[[nn]]
head(dat)
df = merge(df, dat[,c(2:3)],by.x= 1,by.y = 1, all = F)
nam = unique(dat$Set1)
colnames(df)[ncol(df)] = nam
}
# adding HvsM
mol <- read.table("EnrichmentTable_PHG_proteomics_bnGlobalKDA_nodes_sig_over_DEP_HvsM.txt", 
                  sep="\t", stringsAsFactors = FALSE, header = TRUE)
mol<- mol[,c(1,2,9)]
head(mol)
splt <- split(mol, mol$Set1)
for (mm in 1:length(splt)) {
  dat = splt[[mm]]
  head(dat)
  df = merge(df, dat[,c(2:3)],by.x= 1,by.y = 1, all = F)
  nam = unique(dat$Set1)
  colnames(df)[ncol(df)] = nam
}
# adding PCT
mol <- read.table("EnrichmentTable_PHG_proteomics_bnGlobalKDA_nodes_sig_over_PCT.txt", 
                  sep="\t", stringsAsFactors = FALSE, header = TRUE)
mol<- mol[,c(1,2,9)]
head(mol)
splt <- split(mol, mol$Set1)
for (jj in 1:length(splt)) {
  dat = splt[[jj]]
  head(dat)
  df = merge(df, dat[,c(2:3)],by.x= 1,by.y = 1, all = F)
  nam = unique(dat$Set1)
  colnames(df)[ncol(df)] = nam
}

## calculating prot ranks
row.names(df) = df$KDP
mat <- as.matrix(df[,-1])
head(mat)
library(NetWeaver)
source("C:/pc.backup/Documents/R_code/R_functions_for_cor_network/pvalue_aggregation00.R")
# method=c('ACAT','ProductOfRank','MeanOfLog','MeanOfLogLog','Mean')
rel.ACAT = ensemble_rank(mat, method = "ACAT", standardize = T)
rel.ProductOfRank = ensemble_rank(mat, method = "ProductOfRank", standardize = T)
rel.MeanOfLog = ensemble_rank(mat, method = "MeanOfLog", standardize = T)
rel.MeanOfLogLog = ensemble_rank(mat, method = "MeanOfLogLog", standardize = T)
rank = data.frame(ACAT = rel.ACAT, 
                  ProductOfRank = rel.ProductOfRank,
                  MeanOfLog = rel.MeanOfLog,
                  MeanOfLogLog = rel.MeanOfLogLog)
head(rank)
final = merge(df, rank, by.x = 1, by.y = 0, all = F)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed")
load("proteomics.allAdj.RData")
info = protein.info[,c(2,1,6:9)]
head(info)
final$KDP = gsub("\\.Layers3","",final$KDP)
final = merge(info, final, by.x = 1,by.y = 1, all = F)
final <- final[order(final$MeanOfLog, decreasing = T),]
final$rankOrder = 1:dim(final)[1]
row.names(rank) = gsub("\\.Layers3","",row.names(rank))

setwd(wd)
write.table(final, file = "MSBB_PHG_proteomics_globalKDP_BNsubnetwork3layers_rank_04_26_2022.txt",
            sep="\t", row.names = F, quote = F )
write.table(rank, file = "MSBB_PHG_proteomics_globalKDP_BNsubnetwork3layers_rankOnly_04_26_2022.txt",
            sep="\t")
# over
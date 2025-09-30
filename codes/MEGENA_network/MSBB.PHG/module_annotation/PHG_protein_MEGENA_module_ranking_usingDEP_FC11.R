##
rm(list=ls())
wd= "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/hyperGeTest_new"
setwd(wd)

mol <- read.table("EnrichmentTable_Protein_exprCorrected_HvsL_CDR_final_FC11_vs_PHG_protein_MEGENA_modules.txt", 
                      sep="\t", stringsAsFactors = FALSE, header = TRUE)
mol.neg <- mol[mol$Set2 == "negative",c(1,2)]
head(mol.neg)
final.p <- mol.neg[order(mol.neg$Set1),]
replace.df <- final.p 
replace.df$replacement <- rep(1000, dim(replace.df)[1])
replace.df <- replace.df[, c(1,3)]
head(replace.df)


traits = c("HvsL_Braak","HvsL_CDR","HvsL_CERAD","HvsL_PlaqueMean",
           "HvsM_Braak","HvsM_CDR","HvsM_CERAD","HvsM_PlaqueMean",
           "MvsL_CERAD","MvsL_PlaqueMean")

for (nn in 1:length(traits)) {

dat <- read.table(paste0(paste0("EnrichmentTable_Protein_exprCorrected_", traits[nn]),"_final_FC11_vs_PHG_protein_MEGENA_modules.txt"), 
                  sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(dat)
dat00 <- dat[,1:9]
head(dat00)

dat00.neg <- dat00[dat00$Set2 == "negative",c(1,3,9)]
head(dat00.neg)
neg.sort <- dat00.neg[order(dat00.neg$Set1),]
head(neg.sort)
row.names(neg.sort) <- neg.sort$Set1
neg <- neg.sort[,-1]
colnames(neg) <- paste(traits[nn], c("DEG(-).no","FDR(-).p"), sep=".")
head(neg)

dat00.pos <- dat00[dat00$Set2 == "positive",c(1,3,9)]
head(dat00.pos)
pos.sort <- dat00.pos[order(dat00.pos$Set1),]
head(pos.sort)
row.names(pos.sort) <- pos.sort$Set1
pos <- pos.sort[,-1]
colnames(pos) <- paste(traits[nn],c("DEG(+).no","FDR(+).p"), sep=".")
head(pos)

if (dim(neg)[1]== 0) {neg <- replace.df}
if (dim(pos)[1]== 0) {pos <- replace.df}

final <- cbind(neg, pos)
final.p <- cbind(final.p,final[, c(2,4)])
}
head(final.p)
row.names(final.p) <- final.p$Set1
big.table <- final.p[, !grepl("replacement", colnames(final.p))]
head(big.table)
big.table00 <- big.table[,-2]
mat <- as.matrix(big.table00[,-1])
head(mat)

library(NetWeaver)
ProductOfRank <-  ensemble_rank(mat, method= 'ProductOfRank')
MeanOfLog <-  ensemble_rank(mat, method= 'MeanOfLog')
MeanOfLogLog <-  ensemble_rank(mat, method= 'MeanOfLogLog')
Mean <-  ensemble_rank(mat, method= 'Mean')

df <- data.frame(ProductOfRank = ProductOfRank, MeanOfLog = MeanOfLog, MeanOfLogLog = MeanOfLogLog, Mean = Mean)
head(df)
df$Module.id <- row.names(df)
dff <- df[,c(5,1:4)]
head(dff)

dff.mat <- as.matrix(dff[,-1])
cor(dff.mat)

dir = "PHG_protein_mol_ranking"
if(dir.exists(dir)) { print("This dir has already existed!!!")} else {dir.create(dir)}
write.table(big.table00, file="PHG_protein_mol_ranking/PHG_FC11_DEP_vs_MEGENA_moldule_olap_adjP.txt",
            sep="\t", row.names = FALSE, quote = TRUE)
write.table(dff, file="PHG_protein_mol_ranking/PHG_FC11_DEP_based_protein_MEGENA_moldule_ranking.txt",
            sep="\t", row.names = FALSE, quote = TRUE)




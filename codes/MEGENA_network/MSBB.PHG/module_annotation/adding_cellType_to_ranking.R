##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/hyperGeTest_new/PHG_protein_mol_ranking"
setwd(wd)

rank <- read.table("PHG_proteomics_DEP_FC11_based_MEGENA_module_ranking_01_07_2020_updated.txt",
                   sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(rank)
rank00 <- rank

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein/hyperGeTest")

cell <- read.table("EnrichmentTable_MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_cellType.txt",
                   sep="\t", stringsAsFactors = FALSE, header = TRUE)
cell.sel <- cell[, c(1:2,9)]
head(cell.sel)

splt <- split(cell.sel, cell.sel$Set2)
names(splt)

for (nn in 1:length(splt)) {
dat <- splt[[nn]]
dat_ord <- dat[order(dat$Set1), c(1,3)]
row.names(dat_ord) <- dat_ord$Set1
dat_ord <- dat_ord[, -1]
rank00 <- cbind(rank00, dat_ord)
}
colnames(rank00)[(dim(rank00)[2]-5):dim(rank00)[2]] <- paste(c("ast", "end", "mic", "neu", "oli", "opc"), "p.adj", sep="_")

setwd(wd)
write.table(rank00, file = "PHG_proteomics_bigTable.MEGENA_01_07_2020_updated.txt", row.names = FALSE, quote = F, sep="\t")






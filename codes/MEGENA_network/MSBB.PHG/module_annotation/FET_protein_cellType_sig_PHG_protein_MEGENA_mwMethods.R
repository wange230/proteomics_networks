### Enrichment test of gene lists
rm(list=ls())
wd ="C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein"
setwd(wd)

## FET testing olap with MEGENA module

library(class)
library(rpart)

dat <- read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                  sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(dat)
set1 <- dat[,c(2,4)]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

setwd("C:/pc.backup/common_tools_datasets/BrainCellType")
datt <- read.table("top300_human_brainCellType.tsv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(datt)
set2 <- datt[,c(1,2)]
head(set2)
splt2 <- split(set2, set2$cell)
length(splt2)

source("C:/pc.backup/common_tools_datasets/IGAP_signif_p005/Sig_seperate/astrocyte/geneSetOverLap.R")
size = length(union(unique(sort(dat$Symbol)), unique(datt$markers)))  #23002, 39579

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large", BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest"; if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
write.table(enTable, file = "hyperGeTest/EnrichmentTable_MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_cellType.txt",
            sep="\t", row.names = FALSE)


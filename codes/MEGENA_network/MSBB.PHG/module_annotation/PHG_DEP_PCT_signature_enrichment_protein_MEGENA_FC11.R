##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/DEP_PCT_all"

## FET testing olap with MEGENA module
library(class)
library(rpart)

setwd ("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein")
data <- read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                  sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(data)
set1 <- data[,c(1,4)]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

setwd(wd)
dat <- read.table("PCT_FDR_summary.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(dat)
set2 <- dat
head(set2)
splt2 <- split(set2, set2$Module_target)
length(splt2)

source("C:/pc.backup/common_tools_datasets/IGAP_signif_p005/Sig_seperate/astrocyte/geneSetOverLap.R")
size = length(unique(sort(data$Protein)))  #23002, 39579

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest"; if(dir.exists(dir)){print("This dir has already existed!!!")} else{dir.create(dir)}
write.table(enTable, file = "hyperGeTest/EnrichmentTable_all_PCT_vs_PHG_MEGENA_modules.txt", sep="\t", row.names = FALSE)








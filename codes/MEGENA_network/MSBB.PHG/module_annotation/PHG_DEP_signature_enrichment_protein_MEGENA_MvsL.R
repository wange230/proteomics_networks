##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new"
setwd(wd)

## FET testing olap with MEGENA module
library(class)
library(rpart)

setwd(wd)
data <- read.table("MEGENA_reformat_protein/MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                  sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(data)
set1 <- data[,c(1,4)]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019")
dat <- read.table("Protein_exprCorrected_MvsL_Braak_final.tsv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(dat)

dat.sig <- dat[abs(dat$logFC) > log2(1.1) & dat$adj.P.Val < 0.05,]

dat.sig.neg <- dat.sig[dat.sig$logFC < 0,]
dat.sig.pos <- dat.sig[dat.sig$logFC > 0,]

if(dim(dat.sig.neg)[1] == 0 & dim(dat.sig.pos)[1] == 0) {stop("NO DEP existed!!!")}

df <- data.frame(protein.id = c(row.names(dat.sig.neg),row.names(dat.sig.pos)),
                 module = c(rep("negative", dim(dat.sig.neg)[1]),rep("positive", dim(dat.sig.pos)[1]))
                 )
head(df)
set2 <- df
head(set2)
splt2 <- split(set2, set2$module)
length(splt2)

source("C:/pc.backup/common_tools_datasets/IGAP_signif_p005/Sig_seperate/astrocyte/geneSetOverLap.R")
size = length(unique(sort(data$Protein)))  #23002, 39579

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest_new"; if(dir.exists(dir)){print("This dir has already existed!!!")} else{dir.create(dir)}
write.table(enTable, file = "hyperGeTest_new/EnrichmentTable_Protein_exprCorrected_MvsL_Braak_final_FC11_vs_PHG_protein_MEGENA_modules.txt",
            sep="\t", row.names = FALSE)









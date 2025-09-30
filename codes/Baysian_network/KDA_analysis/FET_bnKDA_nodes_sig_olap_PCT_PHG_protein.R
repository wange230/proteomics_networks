##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/MSBB_new_rel"
library(class)
library(rpart)
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/prot_cor_trait2022")
dat <- read.table("PCT2022/phg_prot_cor_traits_signatures.txt",
                  sep="\t",stringsAsFactors = F, header = T)
head(dat)
dat$geneID = paste0(dat$geneID,"_HUMAN")
set1 <- dat[,c(1,2)]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

setwd("C:/Users/wange13/Desktop/MSBB_new_rel/bn_subnet")
datt <- read.table("PHG_protein_KDA_BN_nodes_signatue.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(datt)
set2 <- datt[,c(1,2)]
head(set2)
splt2 <- split(set2, set2$module)
length(splt2)

source("C:/pc.backup/common_tools_datasets/IGAP_signif_p005/Sig_seperate/astrocyte/geneSetOverLap.R")
size = 12147  

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir.create("hyperGeTest2022")
write.table(enTable, file = "hyperGeTest2022/EnrichmentTable_PHG_proteomics_bnGlobalKDA_nodes_sig_over_PCT.txt",
             sep="\t", row.names = FALSE)
## over
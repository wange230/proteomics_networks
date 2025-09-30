##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019/count_DEP_FC11_FDR005_withSymbol"
## FET testing olap with MEGENA module
library(class)
library(rpart)
setwd(wd)
dep <- read.table("DEP_final_02_21_2020_summary.txt", stringsAsFactors = F, sep="\t", header = T)
head(dep)
dat <- dep[,c(2,9)]
set1 <- dat[,c(1,2)]
head(set1)
splt1 <- split(set1, set1$comparisons)
length(splt1)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed")
load("proteomics.allAdj.RData")
gen <- protein.info[,1:2]
head(gen)
comm.sig <- unique(gen$Symbol)

setwd("C:/pc.backup/common_tools_datasets/IGAP_signif_p005")
tbl <- read.table("signature_abeta_adgwas.tsv", sep="\t", stringsAsFactors = F, header = T)
head(tbl)

df <- tbl[tbl$Signature%in%comm.sig,] 
head(df)
datt <- df[,c(2,3)]
set2 <- datt[,c(1,2)]
head(set2)
splt2 <- split(set2, set2$cc)
length(splt2)

source("C:/pc.backup/common_tools_datasets/IGAP_signif_p005/Sig_seperate/astrocyte/geneSetOverLap.R")
size = length(comm.sig) ## 9210

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large", BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest_DEP"; if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
write.table(enTable, file = "hyperGeTest_DEP/EnrichmentTable_DEP_final_FC11_abetaADGWAS.signature.txt",
            sep="\t", row.names = FALSE)




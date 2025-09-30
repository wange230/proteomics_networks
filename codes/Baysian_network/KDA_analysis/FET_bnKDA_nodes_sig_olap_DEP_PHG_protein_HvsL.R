##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/MSBB_new_rel"

traits = c("Braak","CDR","CERAD","PlaqueMean")

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019")
df <- NULL
for (nn in 1:length(traits)) {
DEP <- read.table(paste0(paste0("Protein_exprCorrected_MvsL_", traits[nn]),"_final.tsv"), 
                   sep="\t", stringsAsFactors = FALSE, header = TRUE, quote='"')
DEP <- DEP[DEP$adj.P.Val < 0.05 & abs(DEP$logFC) > log2(1.1),]
DEP$protein.id <- row.names(DEP)
DEP00 <- DEP[DEP$logFC > 0,]
DEP00$trait <- rep(paste(traits[nn],"MvsL_Up",sep="."), dim(DEP00)[1]) 
DEP01 <- DEP[DEP$logFC < 0,]
DEP01$trait <- rep(paste(traits[nn],"_MvsL_Dn",sep="."), dim(DEP01)[1])
DEP <- rbind(DEP00,DEP01)
DEP <- DEP[,c(7,8)]
df <- rbind(df, DEP)
}
head(df)

library(class)
library(rpart)
dat <- df
head(dat)
set1 <- dat[,c(1,2)]
head(set1)
splt1 <- split(set1, set1$trait)
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
write.table(enTable, file = "hyperGeTest2022/EnrichmentTable_PHG_proteomics_bnGlobalKDA_nodes_sig_over_DEP_MvsL.txt",
             sep="\t", row.names = FALSE)
## over
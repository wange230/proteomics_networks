##
rm(list=ls())
wd ="C:/Users/wange13/Desktop/MSBB_new_rel/PHG_proteomics/subnet_new"
setwd(wd)

final <- read.table("MSBB_PHG_proteomics_PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_M5_subnetwork_final00.txt",
                    sep="\t",stringsAsFactors = FALSE, header = TRUE)
head(final)
final$DEP.CDR.color <- ifelse(final$DEP.CDR.sign == "negative", "blue", 
                              ifelse(final$DEP.CDR.sign == "positive", "red",
                                     ifelse(final$DEP.CDR.sign == "neutral", "white", "turquoise")) )
table(final$DEP.CDR.color)

final$DEP.CERAD.color <- ifelse(final$DEP.CERAD.sign == "negative", "blue", 
                              ifelse(final$DEP.CERAD.sign == "positive", "red",
                                     ifelse(final$DEP.CERAD.sign == "neutral", "white", "turquoise")) )
table(final$DEP.CERAD.color)

final$mRNA.CERAD.color <- ifelse(final$mRNA.CERAD.sign == "negative", "blue", 
                                ifelse(final$mRNA.CERAD.sign == "positive", "red",
                                       ifelse(final$mRNA.CERAD.sign == "neutral", "white", "turquoise")) )
table(final$mRNA.CERAD.color)

final$mRNA.CDR.color <- ifelse(final$mRNA.CDR.sign == "negative", "blue", 
                                ifelse(final$mRNA.CDR.sign == "positive", "red",
                                       ifelse(final$mRNA.CDR.sign == "neutral", "white", "turquoise")) )
table(final$mRNA.CDR.color)

names(final)
final00 <- final[, c("Protein","Symbol.2", "Symbol", "is.hub",
                     "DEP.CDR.color","DEP.CERAD.color","mRNA.CDR.color","mRNA.CERAD.color","is.bnGlobalKDA")]
head(final00); final00$Protein.id <- gsub("_HUMAN","",final00$Protein)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/bon_new/Protein_OMS")
meth <- read.table("AD_meth_signatures_proteomics.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
head(meth)
final00$oms.neg <- final00$Protein.id%in%meth[meth$module == "oms.neg",]$protein.id
final00$oms.neg <- ifelse(final00$oms.neg == TRUE, "negative", "none")
final00$oms.pos <- final00$Protein.id%in%meth[meth$module == "oms.pos",]$protein.id
final00$oms.pos <- ifelse(final00$oms.pos == TRUE, "positive", "none")
final00$pos.neg.none <- paste(final00$oms.neg, final00$oms.pos, sep="_")
table(final00$pos.neg.none)
head(final00)

final01 <- final00
final01$DEP.CDR <- rep(1, dim(final01)[1])
final01$DEP.CERAD <- rep(1, dim(final01)[1])
final01$mRNA.CDR <- rep(1, dim(final01)[1])
final01$mRNA.CERAD <- rep(1, dim(final01)[1])

final01$PieChart <- paste0(paste(paste(paste(paste0('colorlist="', final01$DEP.CDR.color), 
                                      final01$DEP.CERAD.color, sep=","), final01$mRNA.CDR.color, sep=","),final01$mRNA.CERAD.color, sep=","), '"')
final01$PieChart[1:2]

final01$PieChart00 <- paste('piechart: attributelist="DEP.CDR,DEP.CERAD,mRNA.CDR,mRNA.CERAD"', final01$PieChart)
final01$PieChart00[1:2]
write.table(final01, file="MEGENA_Network_MSBB_PHG_proteomics_reformat_M5_subnetwork_final00_00.txt", sep="\t", quote = FALSE, row.names = FALSE)



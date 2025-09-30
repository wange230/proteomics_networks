#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/RCG_literature_all_KDP"
setwd(wd)
library(ggplot2)
library(ggpubr)
df = read.table("co_existance_KDP_AD_in_pubmed_withKDPranking.tsv",sep="\t",stringsAsFactors = F, header = T)
head(df)
m3 = read.table("../KDP_in_M3_ranking_literature.txt", sep="\t", stringsAsFactors = F, header = T)
head(m3)

dat = merge(df, m3[,c(1,3,4)], by.x = 1, by.y = 1, all.x = T)
head(dat)
dat[is.na(dat)] = ""

write.table(dat, file="co_existance_KDP_AD_in_pubmed_withKDPranking_plus_M3.tsv", sep="\t",row.names = F, quote = F)
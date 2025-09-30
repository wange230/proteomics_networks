#
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics"
setwd(wd)
library(purrr)
library(dplyr)
library(readxl)
library(ggpubr)
library(readxl)
library(gplots)
library(limma)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(gridExtra)
library(grid)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(stringr)

meta = read.table("MSBB.PHG.proteomics.meta.tsv",
                  sep="\t", stringsAsFactors = F, header = T)
meta$ApoE.1
meta$ApoE.2
meta$ApoE = paste0(meta$ApoE.1,meta$ApoE.2)

info = meta[!is.na(meta$ApoE.1),]

df = read.table("MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj.tsv",
                sep="\t", stringsAsFactors = F, header = T, quote='"')
row.names(df) = df$Protein
datExp = df[,-c(1:9)]
print(table(colnames(datExp) == row.names(meta)))# should all 185 TRUE
prot.meta = df[,c(2,1,6)]
head(prot.meta)

datExp = as.matrix(datExp[,colnames(datExp)%in%row.names(info)])
print(table(colnames(datExp) == row.names(info)))# should all 184 TRUE

# 
table(info$CERAD_final)
table(info$CERAD_final, info$ApoE)

# DEP analysis
library(limma)
#removing 'NA'
info = info[!(is.na(info$CERAD_final)|is.na(info$ApoE)), ]
#info <- info [info$ApoE%in%c("33", "34", "44"),]
#info = info[info$CERAD_final%in%c("NL","defAD"),]
datExp = datExp[, colnames(datExp)%in%row.names(info)]
dim(datExp)
print(table(colnames(datExp) == row.names(info)))# should all 112 TRUE
mat <- as.matrix(datExp)
print(table(colnames(mat) == row.names(info)))
info$APOE_AD = paste(info$ApoE, info$CERAD_final, sep="_")
table(info$APOE_AD)
table(grepl("NL",info$APOE_AD))
table(grepl("defAD",info$APOE_AD)& grepl("34|44", info$APOE_AD))
table(grepl("possAD|probAD",info$APOE_AD)& grepl("34|44", info$APOE_AD))

info$pheno <- ifelse(grepl("NL",info$APOE_AD) ,"NL", 
                     ifelse( grepl("defAD",info$APOE_AD)& grepl("34|44", info$APOE_AD) ,"APOE3444_defAD",
                             ifelse( grepl("possAD|probAD",info$APOE_AD)& grepl("34|44", info$APOE_AD) ,"APOE3444_MCI","others"))) 
info$pheno; table(info$pheno)

#info$pheno <- factor(info$pheno, levels = c("NL", "APOE3444_MCI","APOE3444_defAD","others"))
info$pheno; table(info$pheno)

info00 = info[!info$pheno%in%"others",]
info00$sampleIDs = row.names(info00)

# AHNAK
ahnak = mat[row.names(mat)%in%"sp|Q09666|AHNK_HUMAN",]
table(names(ahnak) == colnames(mat))

df = data.frame(ahnak = as.numeric(ahnak), sampleID = names(ahnak))
head(df)

dff = merge(df, info00, by.x = "sampleID",by.y = "sampleIDs", all = F)
head(dff)
dff$Groups = gsub("def","",dff$pheno)
table(dff$Groups)
dff$Groups = factor(dff$Groups, levels=c("NL","APOE3444_MCI","APOE3444_AD"))
dff$Sex_final = gsub("M","male",gsub("F","female", dff$Sex_final))
dff$gender = str_to_sentence(dff$Sex_final)

#
wdata01 = dff[c("ahnak","Groups","gender")]
head(wdata01)
wdata01$Genotypes = factor(wdata01$Groups, levels = c("NL","APOE3444_MCI","APOE3444_AD"))
head(wdata01)
colnames(wdata01)[1] = "AHANK" 
colnames(wdata01)[3] = "Gender" 
p <- ggboxplot(wdata01, x = "Gender", y = "AHANK",
               color = "Genotypes",
               palette = c("blue","turquoise","red"),size = 1.0, width= 0.6,
               add = "jitter", add.params = list(size = 3.0),
               shape = "Genotypes"
               ) +
  theme(legend.position="none") +
  xlab("") + ylab("AHANK\n (normalized expression)")
p 
pdf(file="DEP_APOE/AHANK_exp_over_APOE_AD_status_by_gender.pdf", height = 3.5, width = 4.5)
print(p)
dev.off()

pp <- ggboxplot(wdata01, x = "Genotypes", y = "AHANK",
               color = "Genotypes",
               palette = c("blue","turquoise","red"),size = 1.0, width= 0.6,
               add = "jitter", add.params = list(size = 3.0),
               shape = "Genotypes"
) +
  theme(legend.position="none",
        axis.title.y =element_text(size=15,face="bold",colour="grey20"),
        axis.text.x = element_text(colour="grey20",
                                   size=12.5,
                                   angle= 30,
                                   hjust= 1.0,
                                   vjust= 1.0,
                                   face="bold"),
        axis.text.y = element_text(colour="grey20",
                                   size=12.5,
                                   #angle= 30,
                                   #hjust= 1.0,
                                   #vjust= 1.0,
                                   face="bold")
        ) +
  xlab("") + ylab("AHANK\n (normalized expression)")
pp 
pdf(file="DEP_APOE/AHANK_exp_over_APOE_AD_status.pdf", height = 4.5, width = 5.0)
print(pp)
dev.off()

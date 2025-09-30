#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/adgwas"
setwd(wd)
library(dplyr)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein")
megena = read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                 sep="\t",stringsAsFactors = F, header = T)

##FET
bkg = unique(megena$Symbol) ## Background size

setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/updated_10172023/v7/excel.files")
library(readxl)
dataset <- read_excel("supplementary.data 4.xlsx", sheet =4, skip=4)
top100 = data_frame(dataset)
top100 = top100[top100$moduleRankingOrder<101,]
top100 = top100[top100$moduleRankingOrder<51,]
# positive set
splt = split(megena, megena$module)
df <- NULL
for(ii in 1:length(splt)){
  tbl = splt[[ii]]%>%distinct(Symbol, .keep_all=T)
  tbl = tbl[!is.na(tbl$Symbol),]
  df = rbind(df, tbl)
}
dat = df[df$Symbol%in%bkg&df$module%in%top100$Module.id,]
#dat = df[df$Symbol%in%bkg,]

set1 <- dat[,c(2,4)]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

## sampling set
# signatures
setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/data_fromMW/ADGWAS_DAA")
daa = read.table("DAA101.txt",sep="\t",stringsAsFactors = F,header = F)
head(daa)
setwd("C:/Users/wange13/Desktop/mouse5fad")
m2h = read.table("Mouse2Human_gene_symbal_conversion_00.txt",
                 sep="\t",stringsAsFactors = F,header = T)
datt = merge(daa, m2h, by.x = 1, by.y =1, all =F)
datt$modules = rep("DAA",nrow(datt))
head(datt)
set2 <- datt[datt$HGNC.symbol%in%bkg,c(2,3)]
head(set2)
splt2 <- split(set2, set2$modules)
length(splt2)

source("c:/pc.backup/Documents/R_code/R_functions_for_cor_network/geneSetOverLap.R")
size = length(bkg)  # 9175

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
write.table(enTable, file = "hyperGeTest/EnrichmentTable_PHG_prot_MEGENA_mod51_for_DAA101.txt",
            sep="\t", row.names = FALSE)
## over
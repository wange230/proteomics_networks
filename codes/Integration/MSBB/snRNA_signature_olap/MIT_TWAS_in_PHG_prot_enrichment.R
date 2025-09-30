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
top100 = data.frame(dataset)
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
adgwas35 = c("ADAR","FRMD4A","SORL1","CCDC50","PSMA1",
             "BIN1","NCK2","PRKD3","GAB2","PTPRJ",
             "SPI1","VASP","RELB","DDB2","NUP160",
             "CKAP5","SCARB2","MS4A6A","APP","PILRA",
             "ABCA1","IGF1R","RIN3","PTK2B","GLIS3",
             "MAN2A1","CASS4","SPON1","TREM2","AMBRA1",
             "MS4A4A","MS4A4E","APOC1","APOE","HLA-DRA")
intersect(adgwas35, bkg)
adgwas35 = data.frame(hgnc_symbol=adgwas35, module="ADGWAS35")

library(readxl)
setwd("C:/Users/wange13/Desktop/new_literature/single-cell/cell_paper")
dataset <- read_excel("mmc1_sun.xlsx",sheet = 14)
adgwas = data.frame(dataset)
adgwas$module = rep("TWAS", nrow(adgwas))
adgwas = adgwas %>% distinct(hgnc_symbol, .keep_all = T)
datt = adgwas[adgwas$hgnc_symbol%in%bkg,c("hgnc_symbol", "module")]
head(datt)
#datt = rbind(adgwas35, datt)
datt$modules = rep("ADWGAS-TWAS", nrow(datt))
set2 <- datt[,c(1,3)]
head(set2)
splt2 <- split(set2, set2$module)
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
write.table(enTable, file = "hyperGeTest/EnrichmentTable_PHG_prot_MEGENA_mod50_for_MIT_TWAS.txt",
            sep="\t", row.names = FALSE)
## over

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
hhm = c("Hexb", "Cst3", "Cx3cr1", "Ctsd",
        "Csf1r", "Ctss", "Sparc", "Tmsb4x", "P2ry12", "C1qa", "C1qb")

dam.1 = c("Tyrobp", "Ctsb", "Apoe", "B2m", "Fth1")

dam.2 = c("Trem2", "Axl", "Cst7", "Ctsl", "Lpl", "Cd9", "Csf1", "Itgax",
          "Clec7a", "Lilrb44", "Timp2")

modc = c("Klk6", "Apod", "Slc5a11", "Pde1a") # mature olig
mfodc = c("Mal", "Mog", "Plp1", "Opalin", "Serinc5", "Ctps1") # mylin forming olig
nfodc = c("Tcf7l2", "Casr", "Cemip2", "Itpr") # newly forming olig

gfap.l.ast = c("Luzp2", "Slc7a10", "Mfge8") #GFAP-low ASC signature
gfap.h.ast = c("Gfap","Id3", "Aqp4", "Myoc", "Id1", "Fabp7")#GFAP-high ASC signature
daa = c("Gfap","Cstb", "Vim", "Osmr", "Gsn", "Ggta1p")

datt = data.frame(genes = c(hhm,dam.1,dam.2,modc,mfodc,nfodc,gfap.l.ast,gfap.h.ast,daa, dam.1,dam.2),
                  modules = c(rep("HHM",length(hhm)),
                              rep("DAM1", length(dam.1)),
                              rep("DAM2", length(dam.2)),
                              rep("OL", length(modc)+length(mfodc)+length(nfodc)),
                              rep("GFAP", length(gfap.l.ast)+length(gfap.h.ast)),
                              rep("DAA", length(daa)),
                              rep("DAM", length(dam.1)+length(dam.2)) )
)
datt$genes = toupper(datt$genes)
head(datt)
set2 <- datt[datt$genes%in%bkg&datt$modules%in%c("DAA","DAM","HHM"),]
set2 <- datt[datt$genes%in%bkg&datt$modules%in%c("GFAP"),]
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
write.table(enTable, file = "hyperGeTest/EnrichmentTable_PHG_prot_MEGENA_mod51_for_GFAP.txt",
            sep="\t", row.names = FALSE)
## over
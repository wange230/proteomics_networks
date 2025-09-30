##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/DEP.validation/"

# determine the background
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics")
df <- read.table("MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj.tsv",
                  sep="\t", stringsAsFactors = F, header = T, quote = '"')
head(df)
msbb.gn = unique(df$Symbol)

setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP")
datt <- read.table("Excluded_NA50_RSOMAP_PFC_proteomics_sex.age_Corrected_HvsL_CERAD.tsv",
                   sep="\t", stringsAsFactors = F, header = T)
head(datt)
ros.gn <- sapply(strsplit(row.names(datt), "\\|"), function(x) x[[1]])
length(unique(ros.gn))

ros.msbb.bkg = intersect(msbb.gn,ros.gn)
length(ros.msbb.bkg) ##7431

# list.files()
tbl <- read.table("MSBB_ROSMAP_protein_signatures_at_protein_GN.level.txt",
                  sep="\t", stringsAsFactors = F, header = T)
head(tbl)

## FET
library(class)
library(rpart)

set1 <- tbl[tbl$cohort == "MSBB",c(1,2)]
head(set1)
set1 = set1[set1$protein.id%in%ros.msbb.bkg,]
splt1 <- split(set1, set1$dn_up)
length(splt1)

set2 <- tbl[tbl$cohort == "ROSMAP",c(1,2)]
head(set2)
set2 = set2[set2$protein.id%in%ros.msbb.bkg,]
splt2 <- split(set2, set2$dn_up)
length(splt2)

source("C:/pc.backup/common_tools_datasets/IGAP_signif_p005/Sig_seperate/astrocyte/geneSetOverLap.R")
size = length(ros.msbb.bkg) #7431

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large", BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest"; 
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}

#write.table(enTable, file = "hyperGeTest/EnrichmentTable_ROSMAP_MSBB_union_unique_DEPs_by_GN.txt",
 #           sep="\t", row.names = FALSE)


## ven diagrams
msbb.dep = tbl[tbl$cohort == "MSBB",c(1,2)]
head(msbb.dep)
msbb.dep.neg <- unique(msbb.dep[msbb.dep$dn_up == "dn",]$protein.id)
msbb.dep.pos <- unique(msbb.dep[msbb.dep$dn_up == "up",]$protein.id)

rosmap.dep = tbl[tbl$cohort == "ROSMAP",c(1,2)]
rosmap.dep.neg <- unique(rosmap.dep[rosmap.dep$dn_up == "dn",]$protein.id)
rosmap.dep.pos <- unique(rosmap.dep[rosmap.dep$dn_up == "up",]$protein.id)

setwd(wd)
dir.create("venn.plots")
#setwd("venn.plots")
wd = paste0(wd, "/","venn.plots")
source("C:/pc.backup/Documents/R_code/R_functions_for_cor_network/Function_4signature_venn_fn_legend_pdf.R")

venn4sig(msbb.dep.neg,msbb.dep.pos,rosmap.dep.neg,rosmap.dep.pos, 
         fn = "MSBBvsROSMAP_DEP","MSBB(Down)","MSBB(Up)","ROSMAP(Down)","ROSMAP(Up)")


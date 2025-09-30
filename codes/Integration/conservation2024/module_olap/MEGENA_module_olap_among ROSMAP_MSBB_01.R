## MSBB prot over MSBB mRNA and ROSMAP mRNA and prot MEGENA module overlap
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA"
setwd(wd)

dir.create("FET.mod_olap")

# *** MSBB prot over MSBB mRNA ***
# msbb prot
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein")
msbb = read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                    sep="\t",stringsAsFactors = F, header = T)
head(msbb)

bkg.msbb = unique(msbb$Symbol)

# msbb mRNA
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/MSBB_RNAseq_latest/MEGENA")
target = read.table("BM_36.multiscale_significant.modules.reformat.tsv",
                    sep="\t",stringsAsFactors = F, header = T)
head(target)
bkg.target = unique(target$Symbol)

# 
bkg = intersect(bkg.target, bkg.msbb)
length(bkg)

# positive set
dat = msbb[msbb$Symbol%in%bkg,]
head(dat)
set1 <- dat[,c(2,4)]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

## sampling set
datt = target[target$Symbol%in%bkg, ]
head(datt)
set2 <- datt[,c(2,4)]
head(set2)
splt2 <- split(set2, set2$Module)
length(splt2)

source("c:/pc.backup/Documents/R_code/R_functions_for_cor_network/geneSetOverLap.R")
size = length(bkg)  # 8367

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
write.table(enTable, file = "FET.mod_olap/EnrichmentTable_MEGENA_modules_MSBB_prot_vs_MSBB_mRNA.txt",
            sep="\t", row.names = FALSE)
# next



rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA"
setwd(wd)

dir.create("FET.mod_olap")
# *** MSBB prot over ROSMAP mRNA ***
# msbb prot
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein")
msbb = read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                  sep="\t",stringsAsFactors = F, header = T)
head(msbb)

bkg.msbb = unique(msbb$Symbol)

# ROSMAP mRNA
setwd("C:/Users/wange13/Desktop/ROSMAP_meth/ROSMAP_MEGENA")
target = read.table("multiscale_significant.modules.reformat.tsv",
                    sep="\t",stringsAsFactors = F, header = T)
head(target)
bkg.target = unique(target$Symbol)

# 
bkg = intersect(bkg.target, bkg.msbb)
length(bkg)

# positive set
dat = msbb[msbb$Symbol%in%bkg,]
head(dat)
set1 <- dat[,c(2,4)]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

## sampling set
datt = target[target$Symbol%in%bkg, ]
head(datt)
set2 <- datt[,c(2,4)]
head(set2)
splt2 <- split(set2, set2$Module)
length(splt2)

source("c:/pc.backup/Documents/R_code/R_functions_for_cor_network/geneSetOverLap.R")
size = length(bkg)  # 8032

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
write.table(enTable, file = "FET.mod_olap/EnrichmentTable_MEGENA_modules_MSBB_prot_vs_ROSMAP_mRNA.txt",
            sep="\t", row.names = FALSE)
# next



## MSBB mRNA and ROSMAP mRNA MEGENA module overlap
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA"
setwd(wd)

dir.create("FET.mod_olap")

# *** MSBB mRNA over ROSMAP mRNA ***
# msbb mRNA
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/MSBB_RNAseq_latest/MEGENA")
msbb = read.table("BM_36.multiscale_significant.modules.reformat.tsv",
                    sep="\t",stringsAsFactors = F, header = T)
head(msbb)
bkg.msbb = unique(msbb$Symbol)

# rosmap mRNA
setwd("C:/Users/wange13/Desktop/ROSMAP_meth/ROSMAP_MEGENA")
target = read.table("multiscale_significant.modules.reformat.tsv",
                    sep="\t",stringsAsFactors = F, header = T)
head(target)
bkg.target = unique(target$Symbol)

# 
bkg = intersect(bkg.target, bkg.msbb)
length(bkg)

# positive set
dat = msbb[msbb$Symbol%in%bkg,]
head(dat)
set1 <- dat[,c(2,4)]
head(set1)
splt1 <- split(set1, set1$Module)
length(splt1)

## sampling set
datt = target[target$Symbol%in%bkg, ]
head(datt)
set2 <- datt[,c(2,4)]
head(set2)
splt2 <- split(set2, set2$Module)
length(splt2)

source("c:/pc.backup/Documents/R_code/R_functions_for_cor_network/geneSetOverLap.R")
size = length(bkg)  # 15413

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
write.table(enTable, file = "FET.mod_olap/EnrichmentTable_MEGENA_modules_MSBB_mRNA_vs_ROSMAP_mRNA.txt",
            sep="\t", row.names = FALSE)
# over
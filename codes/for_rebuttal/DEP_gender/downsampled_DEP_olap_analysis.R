#
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEP_by_sex2024"
setwd(wd)

library(purrr)
library(dplyr)

# female
downsample = read.table("Protein_exprCorrected_defAD_vs_NL_Female_CERAD_final_perm1000.txt",
                        sep="\t", stringsAsFactors = F, header = T)
head(downsample)
table(abs(downsample$logFC)>log2(1.1))
table(downsample$adj.P.Val<0.05)
      
splt = split(downsample, downsample$Protein)
len = map_int(splt, function(x)nrow(x))
# present in > 0.8 of the permuation

splt2 = splt[len>800]
sel = names(splt2)

#up = c()
#dat = splt2[[sel[1]]]
#head(dat)
#sum(dat$logFC>0)
up = map_int(splt2, function(x){sum(x$logFC>0)})
dn = map_int(splt2, function(x){sum(x$logFC<0)})

dep800 = data.frame(Protein = sel, up.num = up, down.num = dn)
head(dep800)
dep800$sign = ifelse(dep800$up.num>800 & dep800$down.num==0, "Up", 
                     ifelse(dep800$down.num>800 & dep800$up.num==0, "Down", "mix"))
table(dep800$sign)

sig1 = dep800$Protein[dep800$sign=="Down"]
sig2 = dep800$Protein[dep800$sign=="Up"]

# male
dep.male = read.table("Protein_exprCorrected_defAD_vs_NL_male_CERAD_final.tsv",
                      sep="\t",stringsAsFactors = F, header = T)
dep.male05 = dep.male[abs(dep.male$logFC)>log2(1.1)&dep.male$adj.P.Val<0.05,]
sig3 = dep.male05$Protein[dep.male05$logFC<0]
sig4 = dep.male05$Protein[dep.male05$logFC>0]

# output sex_specific DEPs
tbl = data.frame(proteinIDs = c(sig1,sig2, sig3,sig4), 
                                "Categories of DEPs" = c(rep("DEPs down regulated in Females", length(sig1)),
                                                       rep("DEPs up regulated in Females", length(sig2)),
                                                       rep("DEPs down regulated in Males", length(sig3)),
                                                       rep("DEPs up regulated in Males", length(sig4)))
                 )

head(tbl)
load("../proteomics.allAdj.RData")

tbl00 = merge(protein.info[,1:2], tbl, by.x = 2, by.y = 1, all = F)
table(tbl00$`Categories of DEPs`)
head(tbl00)
write.table(tbl00, file="Summary_sex_specific_DEPs.txt", sep="\t", row.names = F, quote = F)

# venn 
setwd(wd)
source("C:/pc.backup/Documents/R_code/R_functions_for_cor_network/Function_4signature_venn_fn_legend_pdf.R")

#venn4sig(sig1, sig2, sig3,sig4, "Females.Down","Females.Up", "Males.Down","Males.Up","DEP_female_olap_DEP_male_06_21_2024")

# FET test
# determine background
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed")
load("proteomics.allAdj.RData")
protein = protein.info[,c(1,2,6)]
bkg = protein$Protein ## Background size

# positive set
dat = dep800
head(dat)
set1 <- dat[dat$Protein%in%bkg,c(1,4)]
head(set1)
splt1 <- split(set1, set1$sign)
length(splt1)

## sampling set
datt = dep.male05
head(datt)
datt$module = ifelse(datt$logFC>0, "Up", "Down")
set2 <- datt[datt$Protein%in%bkg,c(1,ncol(datt))]
head(set2)
splt2 <- split(set2, set2$module)
length(splt2)

source("c:/pc.backup/Documents/R_code/R_functions_for_cor_network/geneSetOverLap.R")
size = length(bkg)  # 12147

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
#write.table(enTable, file = "hyperGeTest/EnrichmentTable_DEP_femaleolap_DEP_male.txt",
 #           sep="\t", row.names = FALSE)

## enrichment of olap signatures in MEGENA modules
F_sep_dn = setdiff(sig1, sig3)
F_sep_up = setdiff(sig2, sig4)
F_M_dn = intersect(sig1, sig3)
F_M_up = intersect(sig2, sig4)
M_sep_dn = setdiff(sig3, sig1)
M_sep_up = setdiff(sig4, sig2)

df = data.frame(Protein = c(F_sep_dn,F_sep_up,F_M_dn,F_M_up,M_sep_dn,M_sep_up),
                module = c(rep("Female_sep_dn", length(F_sep_dn)),
                           rep("Female_sep_up", length(F_sep_up)),
                           rep("Female_Male_dn", length(F_M_dn)),
                           rep("Female_Male_up", length(F_M_up)),
                           rep("Male_sep_dn", length(M_sep_dn)),
                           rep("Male_sep_up", length(M_sep_up))
                           
                ) )
head(df)

top30 = c("M3","M182","M475","M750","M34","M39","M193","M492","M35","M210","M245","M762","M46","M510","M771","M533","M9","M221","M520",
          "M5","M70","M769","M43","M306","M659","M806","M2","M120","M319","M186")

# FET test
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein")
megena = read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                    sep="\t",stringsAsFactors = F,header = T)
bkg = unique(megena$Protein) ## Background size

# positive set
dat = df
dat = dat[dat$Protein%in%bkg,]
head(dat)
set1 <- dat[,c(1,ncol(dat))]
head(set1)
splt1 <- split(set1, set1$module)
length(splt1)

## sampling set
datt= megena
datt = datt[datt$Protein%in%bkg,]
head(datt)
set2 <- datt[,c(1,ncol(datt))]
head(set2)
splt2 <- split(set2, set2$module)
length(splt2)

source("c:/pc.backup/Documents/R_code/R_functions_for_cor_network/geneSetOverLap.R")
size = length(bkg)  #  12109

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
dir ="hyperGeTest"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
write.table(enTable, file = "hyperGeTest/EnrichmentTable_PHG_proteomics_modules_olap_spec_DEPs_by_defAD_NL_given_sex_MSBB.txt",
            sep="\t", row.names = FALSE, quote =F)

 

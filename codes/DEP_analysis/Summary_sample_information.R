# summary of sample information
rm(list=ls())

# MSBB
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics")

meta.msbb <- read.table("MSBB.PHG.proteomics.meta.tsv", sep = "\t", stringsAsFactors = F, header = T)
table(meta.msbb$CERAD_final)
table(meta.msbb$CERAD, meta.msbb$CERAD_final)

table(meta.msbb$Sex_final[meta.msbb$CERAD_final == "NL"])
table(meta.msbb$Sex_final[meta.msbb$CERAD_final == "defAD"])
table(meta.msbb$Sex_final[meta.msbb$CERAD_final == "possAD"|meta.msbb$CERAD_final == "probAD"])

summary(meta.msbb$Age_final[meta.msbb$CERAD_final == "NL"])
summary(meta.msbb$Age_final[meta.msbb$CERAD_final == "defAD"])
summary(meta.msbb$Age_final[meta.msbb$CERAD_final == "possAD"|meta.msbb$CERAD_final == "probAD"])
summary(meta.msbb$Age_final)

summary(meta.msbb$PMI_final[meta.msbb$CERAD_final == "NL"])/60
summary(meta.msbb$PMI_final[meta.msbb$CERAD_final == "defAD"])/60
summary(meta.msbb$PMI_final[meta.msbb$CERAD_final == "possAD"|meta.msbb$CERAD_final == "probAD"])/60
summary(meta.msbb$PMI_final)/60

# grouping per trait
table(meta.msbb$Braak_final)
# 0  1  2  3  4  5  6 
# 6 21 24 37 19 19 45 
table(meta.msbb$CDR_final)
#  0   0.5   1   2   3   4   5 
#  28  29   22   26  36  24  20 

table(meta.msbb$PlaqueMean_final==0)
table(meta.msbb$PlaqueMean_final > 9)
table(meta.msbb$PlaqueMean_final > 0 & !meta.msbb$PlaqueMean_final > 9)

#######################
## ROSMAP
rm(list=ls())
setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/output")

meta.rosmap <- read.table("ROSMAP_TNT_proteomics_meta.txt", sep = "\t", stringsAsFactors = F, header = T)
table(meta.rosmap$ceradsc)

table(meta.rosmap$msex[meta.rosmap$ceradsc == 4])
table(meta.rosmap$msex[meta.rosmap$ceradsc == 1])
table(meta.rosmap$msex[meta.rosmap$ceradsc == "3"|meta.rosmap$ceradsc == "2"])

summary(meta.rosmap$age_death[meta.rosmap$ceradsc == "4"])
summary(meta.rosmap$age_death[meta.rosmap$ceradsc == "1"])
summary(meta.rosmap$age_death[meta.rosmap$ceradsc == "3"|meta.rosmap$ceradsc == "2"])
summary(meta.rosmap$age_death)

# grouping per trait
table(meta.rosmap$braaksc)
# 0   1   2   3   4   5   6 
# 6  29  40 120 131  71   3 
table(meta.rosmap$ceradsc)
# 1   2   3   4 
# 114 137  43 106 

# get MMSE 
rm(list=ls())
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/ROSMAP")
meta <- read.csv("AMP-AD_ROSMAP_Rush-Broad_Clinical.csv")
wd="C:/Users/wange13/Desktop/ROSMAP.proteomics/output"
setwd(wd)
dat <- read.table("ROSMAP_TNT_proteomics_expr_sex_aod_adj_lm.txt", sep="\t", stringsAsFactors = F, header = T)
colnames(dat) <- gsub("X","",colnames(dat))
meta<- meta[!is.na(meta$cts_mmse30_lv),]
dat.sel <- dat[,colnames(dat)%in%meta$projid]
meta.sel <- meta[meta$projid%in%colnames(dat.sel),]

mat <- as.matrix(dat.sel[, order(as.numeric(colnames(dat.sel)))])
info <- meta.sel[order(meta.sel$projid),]

table(info$projid == colnames(mat))

info$pheno <- ifelse(info$cts_mmse30_lv > 27, "Low", ifelse(info$cts_mmse30_lv > 17, "Medium","High"))
info$pheno <- factor(info$pheno, levels = c("Low","Medium", "High"))
info$pheno; table(info$pheno)


########################
# BLSA
rm(list=ls())
setwd("C:/Users/wange13/Desktop/Proteomics_lita/jh_proteomics/meta")
meta.blsa = read.csv("blsa_proteomics_pfc_traits.csv")
head(meta.blsa)
meta.blsa47 = meta.blsa[meta.blsa$MFG == 1,]
table(meta.blsa47$MFG)
table(meta.blsa47$Precuneus)
table(meta.blsa47$AGE[meta.blsa47$CT == 1])
table(meta.blsa47$AGE[meta.blsa47$AsymAD == 1])
table(meta.blsa47$AGE[meta.blsa47$AD == 1])

table(meta.blsa47$SEX[meta.blsa47$CT == 1]) # 
#0  1 
#3 10
table(meta.blsa47$SEX[meta.blsa47$AsymAD == 1])
#0  1 
# 4 10 
table(meta.blsa47$SEX[meta.blsa47$AD == 1])
#0  1 
# 10 10 
# over


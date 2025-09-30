##
rm(list=ls())
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/ROSMAP")
meta <- read.csv("AMP-AD_ROSMAP_Rush-Broad_Clinical.csv")

wd="C:/Users/wange13/Desktop/ROSMAP.proteomics"
setwd(wd)
dat <- read.table("output/ROSMAP_TNT_proteomics_expr_sex_aod_adj_lm.txt", sep="\t", stringsAsFactors = F, header = T)
colnames(dat) <- gsub("X","",colnames(dat))
meta<- meta[!is.na(meta$cts_mmse30_lv),]
dat.sel <- dat[,colnames(dat)%in%meta$projid]
meta.sel <- meta[meta$projid%in%colnames(dat.sel),]

mat <- as.matrix(dat.sel[, order(as.numeric(colnames(dat.sel)))])
info <- meta.sel[order(meta.sel$projid),]

table(info$projid == colnames(mat))
# DEP analysis
library(limma)
info$pheno <- ifelse(info$cts_mmse30_lv > 27, "Low", ifelse(info$cts_mmse30_lv > 17, "Medium","High"))
info$pheno <- factor(info$pheno, levels = c("Low","Medium", "High"))
info$pheno; table(info$pheno)

design <- model.matrix(~0+info$pheno)
design
colnames(design) <- c("Low","Medium", "High")
head(design)

matt <- as.matrix(mat) 
matt[1:3,1:3]; dim(matt)
fit <- lmFit(matt, design)
contrast.matrix <- makeContrasts(Medium - Low,  High - Medium, High - Low,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# output
setwd(wd)
outputdir = "DEP.updated"
if(dir.exists("DEP.updated")){print("This dir has already existed!!!")}else{dir.create(outputdir)}

colnames(design)
ml<- topTable(fit2, coef="Medium - Low", number=nrow(fit2))
write.table(ml, file="DEP.updated/Excluded_NA05_RSOMAP_PFC_proteomics_sex.age_corrected_MvsL_MMSE.tsv", sep="\t")
hm<- topTable(fit2, coef="High - Medium", number=nrow(fit2))
write.table(hm, file="DEP.updated/Excluded_NA05_RSOMAP_PFC_proteomics_sex.age_Corrected_HvsM_MMSE.tsv", sep="\t")
hl<- topTable(fit2, coef="High - Low", number=nrow(fit2))
write.table(hl, file="DEP.updated/Excluded_NA05_RSOMAP_PFC_proteomics_sex.age_Corrected_HvsL_MMSE.tsv", sep="\t")
# over

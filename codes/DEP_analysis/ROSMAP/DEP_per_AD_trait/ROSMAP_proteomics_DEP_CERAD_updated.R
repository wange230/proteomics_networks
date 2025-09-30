##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics"
setwd(wd)

meta = read.table("output/ROSMAP_TNT_proteomics_meta.txt",sep="\t",stringsAsFactors = F, header = T) 
meta = meta[!is.na(meta$ceradsc),]
dat = read.table("output/ROSMAP_TNT_proteomics_expr_sex_aod_adj_lm.txt",sep="\t",stringsAsFactors = F, header = T)
colnames(dat) = gsub("X","",colnames(dat))
table(colnames(dat) == meta$projid) # should be 400 TRUE
row.names(meta) = meta$projid

#limfit 
library(limma)
meta$pheno <- ifelse( meta$ceradsc > 3, "Low", 
                               ifelse (meta$ceradsc > 1 , "Medium", "High" ) )
meta$pheno; table(meta$pheno)
meta$pheno <- factor(meta$pheno, levels = c("Low","Medium", "High"))
meta$pheno; table(meta$pheno)

design <- model.matrix(~0+pheno, data = meta)
design
colnames(design) <- c("Low","Medium", "High")
head(design)

matt <- as.matrix(dat) 
matt[1:3,1:3]; dim(matt)
table(colnames(matt) == row.names(design)) # should all 400 TRUE
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
write.table(ml, file="DEP.updated/Excluded_NA05_k50_RSOMAP_PFC_proteomics_sex.age_corrected_MvsL_CERAD.tsv", sep="\t")
hm<- topTable(fit2, coef="High - Medium", number=nrow(fit2))
write.table(hm, file="DEP.updated/Excluded_NA05_k50_RSOMAP_PFC_proteomics_sex.age_Corrected_HvsM_CERAD.tsv", sep="\t")
hl<- topTable(fit2, coef="High - Low", number=nrow(fit2))
write.table(hl, file="DEP.updated/Excluded_NA05_k50_RSOMAP_PFC_proteomics_sex.age_Corrected_HvsL_CERAD.tsv", sep="\t")
# over

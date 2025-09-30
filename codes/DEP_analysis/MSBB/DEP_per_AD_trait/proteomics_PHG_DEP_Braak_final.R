##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed"
setwd(wd)
load("proteomics.allAdj.RData")
protein.expr.syn[1:5,1:5]
protein.info[1:5,1:5]
meta[1:5,1:5]

mat <- as.matrix(protein.expr.syn)
mat[1:5,1:5]

#limfit 
library(limma)

#removing 'NA'
info <- meta[!is.na(meta$Braak_final) == TRUE,]
mat <- mat[,!is.na(meta$Braak_final) == TRUE] 
table(info$SynapseBrainID_final == colnames(mat)) # should all TRUE
info$pheno <- ifelse( info$Braak_final > 4, "High",  ifelse (info$Braak_final > 2 , "Medium", "Low" ) )
info$pheno; table(info$pheno)
info$pheno <- factor(info$pheno, levels = c("Low","Medium", "High"))
info$pheno; table(info$pheno)

design <- model.matrix(~0+info$pheno)
design
colnames(design) <- c("Low","Medium", "High")
head(design)

mat <- as.matrix(mat) 
mat[1:3,1:3]; dim(mat)
fit <- lmFit(mat, design)
contrast.matrix <- makeContrasts(Medium - Low,  High - Medium, High - Low,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# output
setwd(wd)
outputdir = "DEProtein_final_11_12_2019"
if(dir.exists("DEProtein_final_11_12_2019")){print("This dir has already existed!!!")}else{dir.create(outputdir)}

colnames(design)
ml<- topTable(fit2, coef="Medium - Low", number=Inf)
#write.table(ml, file="DEProtein_final_11_12_2019/Protein_exprCorrected_MvsL_Braak_final.tsv", sep="\t")
hm<- topTable(fit2, coef="High - Medium", number=Inf)
#write.table(hm, file="DEProtein_final_11_12_2019/Protein_exprCorrected_HvsM_Braak_final.tsv", sep="\t")
hl<- topTable(fit2, coef="High - Low",  number=Inf)
#write.table(hl, file="DEProtein_final_11_12_2019/Protein_exprCorrected_HvsL_Braak_final.tsv", sep="\t")

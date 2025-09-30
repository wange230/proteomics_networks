##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript"
setwd(wd)
dir.create("results2024")

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed")
load("proteomics.allAdj.RData")
protein.expr.syn[1:5,1:5]
protein.info[1:5,1:5]
meta[1:5,1:5]

mat <- as.matrix(protein.expr.syn)
mat[1:5,1:5]

#limfit 
library(limma)

#removing 'NA'
info <- meta[!is.na(meta$CDR_final) == TRUE,]
mat <- mat[,!is.na(meta$CDR_final) == TRUE] 
table(info$SynapseBrainID_final == colnames(mat)) # should all TRUE
table(info$CDR_final)
info$pheno <- ifelse( info$CDR_final > 0.5, "Demented",  ifelse (info$CDR_final == 0.5, "MCI", "Non_demented" ) )
info$pheno; table(info$pheno)
info$pheno <- factor(info$pheno, levels = c("Non_demented","MCI", "Demented"))
info$pheno; table(info$pheno)

design <- model.matrix(~0+info$pheno)
design
colnames(design) <- c("Non_demented","MCI", "Demented")
head(design)

mat <- as.matrix(mat) 
mat[1:3,1:3]; dim(mat)
fit <- lmFit(mat, design)
contrast.matrix <- makeContrasts(MCI - Non_demented,  Demented - MCI, Demented - Non_demented,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# output
setwd(wd)
#colnames(design)
coefs = colnames(fit2$coefficients)
contrs = gsub(" - ","vs", coefs)

for(ii in 1:length(coefs)){
tbl<- topTable(fit2, coef= coefs[ii], number=Inf)
write.table(tbl, file= paste0("results2024/Protein_exprCorrected_", coefs[ii], ".tsv"),
            sep="\t")
rm(tbl)
}#
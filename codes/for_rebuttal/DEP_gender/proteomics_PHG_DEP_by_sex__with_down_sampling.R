##
rm(list=ls())
wd = "/sc/arion/projects/adineto/permutation_2024/"
setwd(wd)
library(limma)

# input data
df = read.table("MSBB.PHG.proteomics.PMI_Age_race_batch_adj.tsv",
                sep="\t", stringsAsFactors = F, header = T, quote = '"')
row.names(df) = df$Protein
colnames(df)
protein.info = df[,1:9]
mat = df[,10:ncol(df)]
mat <- as.matrix(mat)
mat[1:5,1:5]

meta = read.table("MSBB.PHG.proteomics.meta.tsv",
                  sep="\t", stringsAsFactors = F, header = T)
print(table(row.names(meta) == colnames(mat))) # should be 185 TRUE
colnames(mat) = meta$SynapseBrainID_final

# meta split by sex
meta.f = meta[meta$Sex_final == "F",]
meta.f.nl = meta.f[meta.f$CERAD_final == "NL",]
nrow(meta.f.nl) # 31
meta.f.defad = meta.f[meta.f$CERAD_final == "defAD",]
nrow(meta.f.defad) # 49
meta.m = meta[meta$Sex_final == "M",]
meta.m.nl = meta.m[meta.m$CERAD_final == "NL",]
nrow(meta.m.nl) # 29
meta.m.defad = meta.m[meta.m$CERAD_final == "defAD",]
nrow(meta.m.defad) # 27

# sampling female meta # 1000 times
setwd(wd)
outputdir = "DEP_downsampling"
if(dir.exists("DEP_downsampling")){print("This dir has already existed!!!")
}else{dir.create(outputdir)}

for(ii in 1:1000){
nl.sel <- sample(1:nrow(meta.f.nl), nrow(meta.m.nl))
nl.sel
meta.f.nl.sel = meta.f.nl[nl.sel,]
defad.sel <- sample(1:nrow(meta.f.defad), nrow(meta.m.defad))
defad.sel
meta.f.defad.sel = meta.f.defad[defad.sel,]
info = rbind(meta.f.defad.sel, meta.f.nl.sel)
table(info$CERAD_final)

#limfit 
row.names(info) = info$SynapseBrainID_final
mat00 = mat[,colnames(mat)%in%info$SynapseBrainID_final]
mat00 = mat00[,order(colnames(mat00))]
info = info[order(info$SynapseBrainID_final),]
print(table(info$SynapseBrainID_final == colnames(mat00))) # should all 56 TRUE
info$pheno = info$CERAD_final
info$pheno <- factor(info$pheno, levels = c("NL","defAD"))
info$pheno; table(info$pheno)

design <- model.matrix(~0+pheno, data = info)
design
colnames(design) <- c("NL","defAD")
head(design)

mat00 <- as.matrix(mat00) 
mat00[1:3,1:3]; dim(mat00)
print(table(colnames(mat00) == row.names(design))) # should all 56 TRUE
fit <- lmFit(mat00, design)
contrast.matrix <- makeContrasts(defAD - NL,  levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# output
colnames(design)
coefs = colnames(fit2$contrasts)
contrs = gsub(" - ","_vs_",coefs)
ml<- topTable(fit2, coef = coefs[1], number=Inf)
head(ml)
ml.final = merge(protein.info[,c(1:2)],ml, by.x = 2, by.y = 0, all = F)
ml.final = ml.final[order(ml.final$adj.P.Val),]

write.table(ml.final, file=paste0("DEP_downsampling/Protein_exprCorrected_", contrs[1],"_","Female_CERAD_final_",ii,".tsv"),
             sep="\t", row.names = F, quote = F)
rm(ml.final, ml, mat00,info, fit, fit2, nl.sel,defad.sel)
print(ii)
}# ii
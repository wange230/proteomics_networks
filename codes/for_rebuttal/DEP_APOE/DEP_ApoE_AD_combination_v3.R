## combined cerad into AD (definite, probable, possible), and NL (normal control)
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics"
setwd(wd)
meta = read.table("MSBB.PHG.proteomics.meta.tsv",
                   sep="\t", stringsAsFactors = F, header = T)
meta$ApoE.1
meta$ApoE.2
meta$ApoE = paste0(meta$ApoE.1,meta$ApoE.2)

info = meta[!is.na(meta$ApoE.1),]

df = read.table("MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj.tsv",
                sep="\t", stringsAsFactors = F, header = T, quote='"')
row.names(df) = df$Protein
datExp = df[,-c(1:9)]
print(table(colnames(datExp) == row.names(meta)))# should all 185 TRUE
prot.meta = df[,c(2,1,6)]
head(prot.meta)

datExp = as.matrix(datExp[,colnames(datExp)%in%row.names(info)])
print(table(colnames(datExp) == row.names(info)))# should all 184 TRUE

# 
table(info$CERAD_final)
table(info$CERAD_final, info$ApoE)

# DEP analysis
library(limma)
#removing 'NA'
info = info[!(is.na(info$CERAD_final)|is.na(info$ApoE)), ]
datExp = datExp[, colnames(datExp)%in%row.names(info)]
dim(datExp)
print(table(colnames(datExp) == row.names(info)))# should all 112 TRUE
mat <- as.matrix(datExp)
print(table(colnames(mat) == row.names(info)))
info$APOE_AD = paste(info$ApoE, info$CERAD_final, sep="_")
table(info$APOE_AD)
table(grepl("NL",info$APOE_AD))
table(grepl("defAD",info$APOE_AD)& grepl("34|44", info$APOE_AD))
table(grepl("possAD|probAD",info$APOE_AD)& grepl("34|44", info$APOE_AD))

info$pheno <- ifelse(grepl("NL",info$APOE_AD)& grepl("33", info$APOE_AD),"APOE33_NL", 
                      ifelse( grepl("defAD|possAD|probAD",info$APOE_AD)& grepl("34|44", info$APOE_AD) ,"APOE4_AD",
                                ifelse( grepl("defAD|possAD|probAD",info$APOE_AD) & grepl("33", info$APOE_AD),"APOE33_AD", "others")))
#info$pheno; table(info$pheno)

info$pheno <- factor(info$pheno, 
                     levels = c("APOE33_NL", "APOE33_AD","APOE4_AD","others"))
#info$pheno; table(info$pheno)

design <- model.matrix(~0+info$pheno)
design
colnames(design) <- c("APOE33_NL", "APOE33_AD","APOE4_AD","others")
head(design)

mat <- as.matrix(mat) 
#mat[1:3,1:3]; dim(mat)
table(colnames(mat) == row.names(info))
fit <- lmFit(mat, design)
contrast.matrix <- makeContrasts(APOE4_AD - APOE33_NL,
                                 APOE33_AD - APOE33_NL,
                                 APOE4_AD - APOE33_AD,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# output
setwd(wd)
outputdir = "DEP_APOE2024"
if(dir.exists("DEP_APOE2024")){print("This dir has already existed!!!")}else{dir.create(outputdir)}

colnames(design)
coefs = colnames(fit2$coefficients)
ctrs = gsub(" - ","-vs-",coefs)
for(ii in 1:length(ctrs)){
hl<- topTable(fit2, coef= coefs[ii],  number=Inf)
hl = merge(prot.meta,hl, by.x = 1, by.y = 0, all=F)
write.table(hl, file= paste0("DEP_APOE2024/Protein_exprCorrected_", ctrs[ii],".tsv"),
            sep="\t", row.names = F, quote = F)
rm(hl)
}#ii
## over

# output sample meta data
setwd(wd)
setwd("DEP_APOE2024")
getwd()
#write.table(info, file="PHG_sample_meta._with_APOE_AD.txt", 
 #           sep="\t",quote = F, row.names = F)
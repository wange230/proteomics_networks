##
#remove.packages("Matrix")
#install.packages("Matrix") from tar or zip
#remove.packages("lme4")
#library(Matrix)
#install.packages("lme4", type = "source")
#library(lme4)
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics_lita/jh_proteomics"
setwd(wd)
library(stringr)
library(impute)
library(tidyverse)
library(purrr)
##syn3606086
#list.files()
raw_data = read.csv("AMP-AD_BLSA_EMORY_ProteomicAnalysis_MaxQuantCorrIntensityUniprot.csv")
dim(raw_data)
#colnames(raw_data)
prot.info = raw_data[,1:5]
head(prot.info)
#raw_data$UniqueID
row.names(raw_data) = raw_data$UniqueID
table(raw_data$Reverse) # 0
clean_data = raw_data[is.na(raw_data$Reverse),]
table(raw_data$Potential.contaminant) # 47
clean_data = clean_data[!clean_data$Potential.contaminant == "+",]
table(clean_data$Symbol.or.Gene.Name %in% c(""," ")) # 11
clean_data = clean_data[!clean_data$Symbol.or.Gene.Name %in% c(""," "),]
clean_data = clean_data[, grepl("Intensity",colnames(clean_data))]
#colnames(clean_data)
clean_data = clean_data[,-c(1:2)]
dim(clean_data) ## 4686   47

meta <- read.csv("meta/blsa_proteomics_pfc_traits.csv")
meta$SampleIDs = gsub("_","",paste("Intensity",meta$SanpleName, sep="."))
table(as.character(meta$SampleIDs)%in%colnames(clean_data))

clean_data_pfc = clean_data[, colnames(clean_data)%in%meta$SampleIDs]
meta.pfc = meta[meta$SampleIDs%in%colnames(clean_data_pfc),]
table(colnames(clean_data_pfc) == meta.pfc$SampleIDs)

dir.create("blsa.pfc.final")
save(clean_data_pfc, meta.pfc, file="blsa.pfc.final/PFC47samples_clean_data.RData")

# data normalization
table(meta.pfc$MFG == 1)
table(meta.pfc$Precuneus == 1)

table(meta.pfc$AD)
table(meta.pfc$AsymAD)
table(meta.pfc$CT)
table(meta.pfc$SEX) # 0, 1 3 17, 30
table(meta.pfc$AGE)
table(is.na(meta$PMI))
table(is.na(clean_data_pfc))
# FALSE   TRUE 
# 181692  38550

datExp = as.matrix(clean_data_pfc)
datExp = log2(datExp)
rkeep = rowSums(!is.na(datExp))> 0.5 * ncol(datExp)
table(rkeep) 
#FALSE  TRUE 
#753  3933 
mat = datExp[rkeep,]
mat = as.matrix(mat)
exp <- impute.knn(mat)
dim(exp$data)
exp00 = as.matrix(exp$data)
dim(exp00)
# 3933   47
exp00[1:5,1:5]

#normalization
colMedians <- apply(exp00,2, median) 
exp01 <- sweep(exp00,2,colMedians,'-')
exp01[1:5,1:5]
hist(exp01[,1])
hist(exp01[1,])
table(colnames(exp01) == meta.pfc$SampleIDs)

save(exp01, meta.pfc, file="blsa.pfc.final/PFC47samples_clean_data_normalized.RData")

# covariable correction # age, PMI, sex
load("blsa.pfc.final/PFC47samples_clean_data_normalized.RData")
info = meta.pfc
table(info$SampleIDs == colnames(exp01))
info$SEX = ifelse(info$SEX == 0, "male","female")
info$SEX = factor(info$SEX)
table(info$SEX)
#female   male 
#30     17 
info$AOD = as.numeric(ifelse(info$AGE == ">=90", 90, info$AGE))
#info$AOD.f = ifelse(info$AOD == 90, "high",ifelse(info$AOD>85, "medium","low"))
#info$AOD.f = factor(info$AOD.f)
summary(info$AOD)
summary(info$PMI)

##mixed model
library(lme4)
expr<- apply(exp01,1,function(x, AOD,PMI,SEX){resid(lmer(x ~ AOD+PMI+(1|SEX),REML=F,
                                 control = lme4::lmerControl(calc.derivs = FALSE, 
                                                             check.conv.singular = .makeCC(action = "ignore",  tol = 1e-400),
                                                             check.rankX = "stop.deficient")))+ mean(x)}
             ,PMI = info$PMI, AOD=info$AOD,SEX=info$SEX)
expr_t <- t(expr)
table(colnames(expr_t) == colnames(exp01)) ## check all should TRUE
setwd(wd)
save(expr_t,info, file = "blsa.pfc.final/PFC47samples_clean_data_normalized.RData_LMM_adj.RData")

## DEP2024.final analysis
#limfit 
library(limma)
load("blsa.pfc.final/PFC47samples_clean_data_normalized.RData_LMM_adj.RData")

table(colnames(expr_t) == info$SampleIDs)
info$AD.f <- ifelse(info$CT == 0 & info$AsymAD == 0 & info$AD == 1, "AD",
                       ifelse(info$CT == 0 & info$AsymAD == 1 & info$AD == 0, "MCI",
                              ifelse(info$CT == 1 & info$AsymAD == 0 & info$AD == 0, "NL","none")))
table(info$AD.f)
#AD MCI  NL 
#20  14  13 
info$pheno <- factor(info$AD.f, levels = c("NL","MCI","AD"))
info$pheno; table(info$pheno)

design <- model.matrix(~0+info$pheno)
design
colnames(design) <- c("NL","MCI","AD")
head(design)

mat <- expr_t
mat <- as.matrix(mat) 
table(colnames(mat) == info$SampleIDs)
mat[1:3,1:3]; dim(mat)
fit <- lmFit(mat, design)
contrast.matrix <- makeContrasts(MCI - NL, AD - NL, AD - MCI,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)


# output
setwd(wd)
outputdir = "DEP2024.final"
if(dir.exists("DEP2024.final")){print("This dir has already existed!!!")}else{dir.create(outputdir)}

colnames(design)
ml<- topTable(fit2, coef="MCI - NL", number=nrow(fit2))
head(ml)
ml$protein.id <- row.names(ml)
ml$protein.access <- gsub(".*\\|","",row.names(ml))
ml$GN <- gsub("\\|.*","",row.names(ml))
head(ml)
ml <- ml[,c(7:9,1:6)]
write.table(ml, file="DEP2024.final/PFC47samples_clean_data_normalized_LMM_adj_ASYMvsNL.tsv",
            sep="\t", row.names = F, quote = F)

hm<- topTable(fit2, coef="AD - MCI", number=nrow(fit2))
hm$protein.id <- row.names(hm)
hm$protein.access <- gsub(".*\\|","",row.names(hm))
hm$GN <- gsub("\\|.*","",row.names(hm))
hm <- hm[,c(7:9,1:6)]
write.table(hm, file="DEP2024.final/PFC47samples_clean_data_normalized_LMM_adj_ADvsASYM.tsv",
            sep="\t", row.names = F, quote = F)

hl<- topTable(fit2, coef="AD - NL", number=nrow(fit2))
hl$protein.id <- row.names(hl)
hl$protein.access <- gsub(".*\\|","",row.names(hl))
hl$GN <- gsub("\\|.*","",row.names(hl))
hl <- hl[,c(7:9,1:6)]
write.table(hl, file="DEP2024.final/PFC47samples_clean_data_normalized_LMM_adj_ADvsNL.tsv",
            sep="\t", row.names = F, quote = F)

# check the results
head(hl)
hl05 = hl[hl$adj.P.Val<0.05,]
length(unique(hl05$GN))
length(unique(hl05$GN[hl05$logFC>0]))
length(unique(hl05$GN[hl05$logFC<0]))

lt = list(ml, hm, hl)
ctrs = c("ASYMvsNL","ADvsASYM","ADvsNL")
for(jj in 1:length(lt) ){
tbl = lt[[jj]]
splt = split(tbl, tbl$GN)
df = NULL
for(ii in 1:length(splt)){
datt = splt[[ii]] 
if(nrow(datt)>1){datt = datt[order(datt$P.Value, decreasing = F),]; df = rbind(df, datt[1,])} else {
  df = rbind(df, datt)
}
}
head(df)
write.table(df, file= paste0("DEP2024.final/PFC47samples_clean_data_normalized_LMM_adj_", ctrs[jj],"_unique_prot.tsv"),
            sep="\t", row.names = F, quote = F)
rm(df)
rm(splt, tbl)
}
# over

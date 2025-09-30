##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics"
setwd(wd)

gen <- read.csv("input/Uniprot2Symbol_Human_v03_SimpleforR.csv")
head(gen)

meta <- read.csv("input/rosmap_50batch_specimen_metadata_for_batch_correction.csv")
head(meta)
meta.sel <- meta[!is.na(meta$projid),]

setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed")
dat <- read.table("C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400_50na.txt",
                  sep="\t", stringsAsFactors = F, header = T)
dat[1:5,1:5]
dattt <- dat

df <- dattt[, colnames(dattt)%in%as.character(meta.sel$SampleID)]
df.ord <- df[, order(colnames(df))]
meta.ord <- meta.sel[order(meta.sel$SampleID),]
table(colnames(df.ord) == meta.ord$SampleID)  ## should all TRUE
colnames(df.ord) <- meta.ord$projid

#ROSMAP-RNAseq analysis
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/ROSMAP/RNAseq/covariates_adjustment")
meta00 <-read.csv("meta.ROSMAP.subtypes.matchedMSBB.121919.tsv", header = TRUE, sep ='\t')
meta00 <- meta00[,c(3,8,17)]
setwd("C:/pc.backup/AD_grants/misc/Data_from_RUSH/ROSMAP_sample_selection")
info <- read.table("dataset_534_basic.tsv", sep="\t", stringsAsFactors = F, header = T)
#info00 <- merge(meta00, info, by.x = 1, by.y = 1, all = FALSE)
info00 = info
dff0 <- df.ord[, colnames(df.ord)%in%info00$projid]
dff0_ord <- dff0[,order(as.numeric(colnames(dff0)))]
metaa.sel <- info00[info00$projid%in%colnames(df.ord),]
metaa.sel.ord <- metaa.sel[order(metaa.sel$projid),]
table(colnames(dff0_ord) == metaa.sel.ord$projid)
#setwd(wd)
# write.table(metaa.sel.ord, file="ROSMAP_TNT_proteomics_meta.txt", sep="\t", quote = F, row.names = F)

metaa.sel.ord$sex.f <- factor(metaa.sel.ord$msex)
metaa.sel.ord$aod <- metaa.sel.ord$age_death
mat <- na.omit(as.matrix(dff0_ord))
expr <- apply(mat, 1, function(x, aod, sex.f){resid(lm(x ~ aod+sex.f))+mean(x)},
              sex.f = metaa.sel.ord$sex.f, aod = metaa.sel.ord$aod)
dim(expr)
expr_t <- t(expr)
dim(expr_t) # protein id at row.names
setwd(wd)
# write.table(expr_t, file="ROSMAP_TNT_proteomics_expr_sex_aod_adj_lm.txt", sep="\t", quote = F)

#limfit 
library(limma)
metaa.sel.ord$pheno <- ifelse( metaa.sel.ord$braaksc > 4, "High", 
                               ifelse (metaa.sel.ord$braaksc > 2 , "Medium", "Low" ) )
metaa.sel.ord$pheno; table(metaa.sel.ord$pheno)
metaa.sel.ord$pheno <- factor(metaa.sel.ord$pheno, levels = c("Low","Medium", "High"))
metaa.sel.ord$pheno; table(metaa.sel.ord$pheno)

design <- model.matrix(~0+metaa.sel.ord$pheno)
design
colnames(design) <- c("Low","Medium", "High")
head(design)

matt <- as.matrix(expr_t) 
matt[1:3,1:3]; dim(matt)
fit <- lmFit(matt, design)
contrast.matrix <- makeContrasts(Medium - Low,  High - Medium, High - Low,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# output
setwd(wd)
outputdir = "DEP.updated"
if(dir.exists("DEP.updated")){print("This dir has already existed!!!")}else{dir.create(outputdir)}

#colnames(design)
ml<- topTable(fit2, coef="Medium - Low", number=nrow(fit2))
write.table(ml, file="DEP.updated/Excluded_NA50_RSOMAP_PFC_proteomics_sex.age_corrected_MvsL_Braak.tsv", sep="\t")
hm<- topTable(fit2, coef="High - Medium", number=nrow(fit2))
write.table(hm, file="DEP.updated/Excluded_NA50_RSOMAP_PFC_proteomics_sex.age_Corrected_HvsM_Braak.tsv", sep="\t")
hl<- topTable(fit2, coef="High - Low", number=nrow(fit2))
write.table(hl, file="DEP.updated/Excluded_NA50_RSOMAP_PFC_proteomics_sex.age_Corrected_HvsL_Braak.tsv", sep="\t")
# over
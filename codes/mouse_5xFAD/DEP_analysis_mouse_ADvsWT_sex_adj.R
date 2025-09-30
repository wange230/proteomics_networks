##
rm(list=ls())
wd="C:/Users/wange13/Desktop/mouse5fad"
setwd(wd)

meta <- read.table("meta_5fad_protein.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(meta)
meta$Genotype <- gsub("5x","",meta$Genotype)
meta$sex.f <- factor(meta$Sex)
sel = c(3,6,12)
df <- read.csv("5xFAD_mouse_wholeProteome_PengLab_at_StJude_v1.1.0.csv")
row.names(df) <- df$Protein.Accession
gen <- df[,1:2]; head(gen)
df.sel <- df[, c(4:dim(df)[2])]
df.completed <- na.omit(df.sel)
table(meta$Sample.ID == colnames(df.completed)) ## should all TRUE
matt <- as.matrix(df.completed)
mat.log2 <- log2(matt)
expr <- apply(mat.log2, 1, function(x, sex.f){resid(lm(x ~ sex.f))+mean(x)},sex.f = meta$sex.f)
dim(expr)
expr_t <- t(expr)
dim(expr_t)

# write.table(expr_t, file ="5xFAD_mouse_wholeProteome_PengLab_at_StJude_v1.1.0.sex.adj.txt", sep="\t",quote = F)

for (nn in 1:length(sel)) {
meta.3 <- meta[meta$Age..mo. == sel[nn],]
head(meta.3)
mat <- as.matrix(expr_t[, colnames(expr_t)%in% meta.3$Sample.ID])

## DEP analysis ##limma analysis
library(limma)
pheno <- factor(meta.3$Genotype, levels = c("WT", "FAD")); length(pheno)
design <- model.matrix(~0+pheno)
head(design)
colnames(design) <- c("WT", "FAD")
head(design)

fit <- lmFit(mat, design)
contrast.matrix <- makeContrasts(FAD - WT,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = TRUE)

setwd(wd)
dir = "mouse_FAD_vs_WT_DEP"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}

colnames(design)
ml <- topTable(fit2, coef="FAD - WT", number=nrow(fit2))
head(ml)
final <- merge(gen, ml, by.x = 2, by.y = 0,all = FALSE)
head(final)

write.table(final, file= paste0(paste0("mouse_FAD_vs_WT_DEP/mouse_FAD_vs_WT_month", sel[nn]), "_Sex_adj_DEP.tsv"),
            sep="\t", row.names = FALSE)
}

##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed"
setwd(wd)
load("proteomics.allAdj.RData")

table(colnames(protein.expr) == row.names(meta))
table(colnames(protein.expr.syn) == meta$SynapseBrainID_final)
names(meta)
protein.expr.complete <- na.omit(protein.expr)
protein.expr.complete[1:5,1:5]

traits <- meta[, c("Braak_final","CDR_final","CERAD_final","PlaqueMean_final")]
dim(traits); table(traits$CERAD_final)
traits$CERAD_final00 <- ifelse(traits$CERAD_final == "NL",1, ifelse(traits$CERAD_final == "possAD", 2, 
                                 ifelse(traits$CERAD_final == "probAD", 3,4)))
table(traits$CERAD_final00)
head(traits)

traits.t <- t(traits[,c(1,2,5,4)])
traits.t[,1:5]

table(colnames(traits.t) == colnames(protein.expr.complete))

mat <- as.matrix(rbind(traits.t, protein.expr.complete))
table(mat[1,])

targets <- data.frame(traits = c("Braak_final","CDR_final","CERAD_final00","PlaqueMean_final"), category = c(1,1,1,0))
head(targets)


library(groupCor)
#groupCor::one2All_cor
#groupCor:::all2all.cor

rels.cor <- groupCor::all2all.cor(mat, targets)
rels.cor$geneID <- paste0(rels.cor$geneID,"_HUMAN")
head(rels.cor)
pro <- protein.info[,c(1:2)]
head(pro)
rels.cor.id <- merge(pro, rels.cor, by.x = 2, by.y = 1, all.y = TRUE)
head(rels.cor.id)

setwd(wd)
dir = "PHG_traits_cor"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}

#write.table(rels.cor, file = "PHG_traits_cor/proteomics.allAdj_cor_traits_PHG.txt",
 #           sep="\t", row.names = FALSE, quote = F)
#write.table(rels.cor.id, file = "PHG_traits_cor/proteomics.allAdj_cor_traits_PHG_withID.txt",
 #           sep="\t", row.names = FALSE, quote = F)

rels.cor[!grepl("_HUMAN", rels.cor$geneID),]
head(rels.cor)







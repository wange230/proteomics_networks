##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics"
setwd(wd)
protein <- read.table("MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj.tsv", 
                      sep="\t", stringsAsFactors = FALSE, header = TRUE, quote = '"')
names(protein)
length(unique(protein$Symbol));length(protein$Symbol)
length(unique(protein$HGNC.EnsemblID));length(protein$HGNC.EnsemblID)
length(unique(protein$Protein));length(protein$Protein)
protein[2,1:9]
row.names(protein) <- protein$Protein
protein.info <- protein[,1:9]
protein.expr <- protein[,10:(dim(protein)[2])]

meta <- read.table("MSBB.PHG.proteomics.meta.tsv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
table(gsub("\\.","-",row.names(meta)) == meta$RunName)
gsub("\\.","-",row.names(meta))[!gsub("\\.","-",row.names(meta)) == meta$RunName]
meta$RunName[!gsub("\\.","-",row.names(meta)) == meta$RunName]
table(meta$SynapseBrainID == meta$SynapseBrainID_final)

table(row.names(meta) == colnames(protein.expr))
protein.expr.syn = protein.expr
table(row.names(meta) == colnames(protein.expr.syn))
colnames(protein.expr.syn) <- meta$SynapseBrainID_final

dir = "processed"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
save(protein.expr, protein.expr.syn,protein.info, meta, file= "processed/proteomics.allAdj.RData")




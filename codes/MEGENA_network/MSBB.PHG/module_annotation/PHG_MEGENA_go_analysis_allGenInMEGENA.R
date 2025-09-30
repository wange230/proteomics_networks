##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein"
setwd(wd)
dat <- read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                   sep="\t", stringsAsFactors = F, header = T)
head(dat)
setwd("C:/pc.backup/common_tools_datasets/Gene_ontology/DB/")
go <- read.table("MSigDB_Selection_v6_1.background_gene_list.tsv", sep="\t", stringsAsFactors = F, header = T)
head(go)
dat.sel <- dat[dat$Symbol%in%go$gene.name,]
gen = unique(sort(dat.sel$Symbol))
df <- data.frame(GeneSymbol = gen, module = rep("combined", length(gen)) )
head(df)
df00 <- data.frame(GeneSymbol = gen, module = rep("combined00", length(gen)) )
head(df00)
dff <- rbind(df,df00)
setwd(wd)
write.table(dff, file = "MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_GenSelModuleAll.tsv",
            sep="\t", row.names = F, quote = F)
## Genontology test for astrocytes_DEG
# enrichment test of gene lists
#*************************************

library(class)
library(rpart)

setwd(wd)
source("C:/pc.backup/common_tools_datasets/Gene_ontology/R_go_functions.R")

# minputfnames = "mouse_inflammatome_signatures.xls"
# idxGS=3; # which column holds gene symbols
# ontologyfname = "DB/GO-MouseV44k-16555gsymbols_trimed.xls" # mouse genome annotation
# enrichType= "GO"

minputfnames = "MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_GenSelModuleAll.tsv";
idxGS= 1; # which column holds gene symbols ##original idxGS=1
#ontologyfname = "Gene_ontology/DB/GO-HumanV3-16773gsymbols_trimed.xls" # human genome annotation
#enrichType= "GO"

ontologyfname = "C:/pc.backup/common_tools_datasets/Gene_ontology/DB/MSigDB.Selection.v6.1.tsv"
enrichType= "MSigDB_updated"

moutdir = paste(enrichType, "/", sep="")
dir.create(moutdir)

moduleBasedIntersectionMatrix_GeneOntology(fnames=minputfnames, fnameB=ontologyfname, outputDir=moutdir, 
                                           uniqueIdCol=idxGS, ctopNrows=10, signifLevel=0.005, removeDuplicate=TRUE, removeGrey=TRUE)
############################## END #################################################

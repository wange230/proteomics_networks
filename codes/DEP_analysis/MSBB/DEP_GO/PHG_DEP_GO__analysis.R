##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019/count_DEP_FC11_FDR005_withSymbol"
setwd(wd)

dat <- read.table("DEP_final_02_21_2020_summary.txt", sep="\t", stringsAsFactors = F, header = T)
head(dat)
setwd("C:/pc.backup/common_tools_datasets/Gene_ontology/DB/")
go <- read.table("MSigDB_Selection_v6_1.background_gene_list.tsv", sep="\t", stringsAsFactors = F, header = T)
head(go)
dat.sel <- dat[dat$Symbol%in%go$gene.name,]

setwd(wd)
write.table(dat.sel, file = "DEP_final_02_21_2020_summary_for_go.txt", sep="\t", row.names = F, quote = F)

## Genontology test for astrocytes_DEG
# enrichment test of gene lists
#*************************************

library(class)
library(rpart)

setwd(wd)
source("Gene_ontology/R_go_functions.R")

# minputfnames = "mouse_inflammatome_signatures.xls"
# idxGS=3; # which column holds gene symbols
# ontologyfname = "DB/GO-MouseV44k-16555gsymbols_trimed.xls" # mouse genome annotation
# enrichType= "GO"

minputfnames = "DEP_final_02_21_2020_summary_for_go.txt";
idxGS= 2; # which column holds gene symbols ##original idxGS=1
#ontologyfname = "Gene_ontology/DB/GO-HumanV3-16773gsymbols_trimed.xls" # human genome annotation
#enrichType= "GO"


ontologyfname = "Gene_ontology/DB/MSigDB.Selection.v6.1.tsv"
enrichType= "MSigDB"

moutdir = paste(enrichType, "/", sep="")
dir.create(moutdir)

moduleBasedIntersectionMatrix_GeneOntology(fnames=minputfnames, fnameB=ontologyfname, outputDir=moutdir, 
                                           uniqueIdCol=idxGS, ctopNrows=10, signifLevel=0.005, removeDuplicate=TRUE, removeGrey=TRUE)
############################## END #################################################

##
rm(list=ls())
setwd("C:/Users/wange13/Desktop/mouse5fad")
source("C:/pc.backup/common_tools_datasets/Mou2Hu_gene_conversion_function.R")

dat <- read.csv("5xFAD_mouse_transcriptome_PengLab_at_StJude_v1.1.0.csv")
musGenes=as.character(dat$GN)
mo2hu.homologue <- convertMouseGeneList(musGenes)
dim(mo2hu.homologue); head(mo2hu.homologue)

write.table(mo2hu.homologue, file= "Mouse2Human_gene_symbal_conversion_00.txt",
            sep="\t", row.names = FALSE)
#convertHumanGeneList("FGF1")

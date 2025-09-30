#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript"
setwd(wd)
dir.create("signal_maps")
library(dplyr)
library(purrr)
library(clusterProfiler)
library(DOSE)
library(enrichplot)

# input go terms
setwd("data_fromMW/go_signal_map")
grp = readRDS("go.hierarchy.RDS")
head(grp)

go = readRDS("go.terms.RDS")
head(go)
go.bp = go[go$Ontology == "BP",]
term2gene = go[,c(2,3)]
head(term2gene)

# DEP signatures
setwd(wd)
dep = read.table("../ahnak_proteomics/GSEA/AHNAK_DEP07_for_GSEA.txt",
                 sep="\t",stringsAsFactors = F, header = T)
dep = dep[dep$P.Value<0.30,]
dep = dep[order(dep$logFC, decreasing = T),]
head(dep)
genelist = dep$logFC
names(genelist) = dep$GN
genelist

gsea = GSEA(
            geneList = genelist,
            exponent = 1,
            minGSSize = 10,
            maxGSSize = 500,
            eps = 1e-30,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            TERM2GENE = term2gene,
            verbose = TRUE,
            seed = 12345678, # FALSE,
            by = "fgsea"
            )
head(gsea)

rels = gsea@result
head(rels)

setwd(wd)
save(rels, file = "signal_maps/AHNAK_DEP03_for_GSEA_BPonly_11012023.RData")
# over
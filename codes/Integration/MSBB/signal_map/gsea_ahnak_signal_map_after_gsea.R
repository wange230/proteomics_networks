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
term2gene = go[,c(2,3)]
head(term2gene)

table(go$go_id %in% union(grp$parent,grp$child))
table(term2gene$go_id %in% union(grp$parent,grp$child))
#FALSE    TRUE 
#125035 4119785

setwd(wd)
list.files("signal_maps")
load("signal_maps/updated/AHNAK_DEP03_for_GSEA_BPonly_11012023.RData")
# 
rels01 = rels[rels$qvalue<0.05,]
colnames(rels01)
table(rels01$NES>0)
head(rels01)
colnames(rels01)
rels01 = rels01[c("ID","Description","NES","qvalue")]
head(go)
go.terms = go %>% distinct(go_id, .keep_all = T)
go.terms = go.terms[,c(2,4,5)]
head(go.terms)
go.terms = go.terms[go.terms$Ontology == "BP",]

# project the nodes in rels01 to the GO net
head(rels01)
head(grp)
# adding parent 
df = merge(grp, rels01, by.x = 1, by.y = 1, all = F )
head(df)
# adding child
dff = merge(df, rels01, by.x = 2, by.y = 1, all = F)
colnames(dff) = gsub("\\.y",".child",colnames(dff)) 
colnames(dff) = gsub("\\.x",".parent",colnames(dff)) 
head(dff)
table(dff$child == dff$Description.child) # should all TRUE
table(dff$parent == dff$Description.parent) # should all TRUE

length(union(dff$child,dff$parent))
length(unique(union(dff$child,dff$parent)))
# 
setwd(wd)
write.table(dff, file = "signal_maps/AHNAK_DEP03_for_GSEA_BP_11012023.txt",
            sep="\t", row.names = F, quote = F)
# over
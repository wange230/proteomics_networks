#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/signal_maps/updated"
setwd(wd)
library(dplyr)
library(purrr)

go.terms = readRDS("../../data_fromMW/go_signal_map/go.terms.RDS")
head(go.terms)
go.term = go.terms %>% distinct(go_id, .keep_all = T)
head(go.term)
go.bp = go.term[go.term$Ontology == "BP",c(2,4,5)]
head(go.bp)

df = read.table("AHNAK_DEP03_for_GSEA_BP_11012023.txt",
                 sep="\t", stringsAsFactors = F,header = T)
head(df)
# remove nodes having no term description/annotation
dff = df[df$child%in%go.bp$go_id,]
dfff = dff[dff$parent%in%go.bp$go_id,]
write.table(dfff, file="AHNAK_DEP03_for_GSEA_BP_removing_go_no_annot_11012023.txt",
            sep="\t", row.names = F, quote = F)

# get the nodes
nodes = data.frame(nodes = union(dff$parent, dff$child) )

head(nodes)
nodes$isParent = ifelse(nodes$nodes %in% df$parent, "Yes", "No")
nodes$isChild  = ifelse(nodes$nodes %in% df$child, "Yes", "No")
head(nodes)

nodes$parentOnly = ifelse(nodes$isParent == "Yes" & nodes$isChild == "No", "Yes","No")
nodes$childOnly = ifelse(nodes$isParent == "No" & nodes$isChild == "Yes", "Yes","No")
nodes$both = ifelse(nodes$isParent == "Yes" & nodes$isChild == "Yes", "Yes","No")



nodes00 = merge(nodes, go.bp, by.x = 1, by.y = 1, all = F)
head(nodes00)
table(nodes00$both)

load("AHNAK_DEP03_for_GSEA_BPonly_11012023.RData")
# 
rels01 = rels[rels$qvalue<0.05,]
table(rels01$NES>0)
head(rels01)
rels01 = rels01[c("ID","Description","NES","qvalue")]

nodes01 = merge(nodes00, rels01, by.x = 1, by.y = 1, all = F)
nodes01$FDR = log10(nodes01$qvalue) * (-1)
summary(nodes01$FDR)
nodes01$FDR.cat = ifelse(nodes01$FDR >summary(nodes01$FDR)[5], "High",
                         ifelse(nodes01$FDR >summary(nodes01$FDR)[2], "Medium", "Low"))
table(nodes01$FDR.cat)
nodes01$NES.sign = ifelse(nodes01$NES>0, "Activated", 
                          ifelse(nodes01$NES <0, "Suppressed", "none"))

write.table(nodes01, file="AHNAK_DEP03_for_GSEA_BP_nodes_11012023.txt",
             sep="\t", row.names = F, quote = F)
# over
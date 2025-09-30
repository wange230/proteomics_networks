##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/MSBB_new_rel/PHG_proteomics/subnet_new"
setwd(wd)
megena <- read.table("MSBB_PHG_proteomics_PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_M3_subnetwork.txt",
                      sep="\t", stringsAsFactors = F, header = T)
head(megena)
megena$Source <- megena$col
megena$Target <- megena$row
megena$Weight <- rep(1,dim(megena)[1])
megena$Type <- rep("Undirected",dim(megena)[1])
head(megena)

notes <- union(megena$Source,megena$Target)
nodes <- data.frame(Id = 1:length(notes), Label = notes)
head(nodes)

setwd(wd)
dep <- read.table("MSBB_PHG_proteomics_PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_M3_subnetwork_nodes.txt",
                   sep="\t", stringsAsFactors = F, header = T)
head(dep)
dep = dep[,c(2:3,6:10)]
head(dep)

nodes00 <- merge(nodes, dep, by.x = 2, by.y = 1, all.x = T)
head(nodes00)
nodes.final <- nodes00[,c(2,1,3:8)]
nodes.final00 <- nodes.final
nodes.final00$Label = ifelse(!nodes.final00$is.bnGlobalKDA == "no", nodes.final00$Symbol, "")


megena.source <- merge(megena[,3:6], nodes, by.x = 1, by.y = 2, all.x = T)
head(megena.source)
megena.source00 <- megena.source[,c(5,2:4)]
colnames(megena.source00)[1] <- "Source"
head(megena.source00)

megena.target <- merge(megena.source00, nodes, by.x = 2, by.y = 2, all.x = T)
head(megena.target)
megena.target00 <- megena.target[,c(2,5,4,3)]
colnames(megena.target00)[2] <- "Target"
megena.gephi <- megena.target00[order(megena.target00$Source),]
head(megena.gephi)

setwd(wd)
dir = "gephi.data"; dir.create(dir)
write.csv(nodes.final00, file = "gephi.data/global.M3subMEGENA_nodes_gephi.csv", quote = F, row.names = F)
write.csv(megena.gephi, file = "gephi.data/global.M3subMEGENA_edges_gephi.csv", quote = F, row.names = F)


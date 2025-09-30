##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new"
setwd(wd)
megena <- read.table("MEGENA_Network.txt", sep="\t", stringsAsFactors = F, header = T)
megena$Source <- megena$row
megena$Target <- megena$col
megena <- megena[,c(4,5)]
head(megena)
megena$Weight <- rep(1,dim(megena)[1])
megena$Type <- rep("Undirected",dim(megena)[1])
head(megena)

notes <- union(megena$Source,megena$Target)
nodes <- data.frame(Id = 1:length(notes), Label = notes)
head(nodes)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/MEGENA_reformat_protein")
dep <- read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv",
                   sep="\t", stringsAsFactors = F, header = T)
mol.sel = c("M2","M3","M4","M5","M6","M7","M8","M9","M10","M11","M12","M13","M15","M16","M17","M18","M19","M20")
#dep = dep[dep$module%in%c("M3","M5","M9"),]
dep = dep[dep$module%in%mol.sel,]
head(dep)

nodes00 <- merge(nodes, dep, by.x = 2, by.y = 1, all.x = T)
head(nodes00)
nodes.final <- nodes00[,c(2,1,4:5)]
head(nodes00)
nodes00$is.hub <- ifelse(is.na(nodes00$is.hub), "",nodes00$is.hub)
nodes00$module <- ifelse(is.na(nodes00$module), "",nodes00$module)
head(nodes00)

megena.source <- merge(megena, nodes, by.x = 1, by.y = 2, all.x = T)
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
write.csv(nodes.final, file = "gephi.data/global.megena_nodes_gephi_1stlayers.csv", quote = F, row.names = F)
write.csv(megena.gephi, file = "gephi.data/global.megena_edges_gephi_1stlayers.csv", quote = F, row.names = F)


### bn networks finding downstream targets
rm(list=ls())
setwd("C:/Users/wange13/Desktop/MSBB_new_rel/hyperGeTest/PHG_protein_bnGlobalKDA_ranking")
kda <- read.table("PHG_FC11_DEP_based_protein_PHG_protein_bnGlobalKDA_rankingWithID.txt",
                  sep="\t", stringsAsFactors = FALSE, header = TRUE)
kda <- kda[order(kda$MeanOfLog, decreasing = T),]
head(kda)

wd <- "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/bn_network"
setwd(wd)
library(KDA)
ntwk <- read.table("global.BN_cys.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)
head(ntwk)
colnames(ntwk) <- c("from","to")
ntwk$from <- gsub("_HUMAN","",gsub(".*:","",ntwk$from))
ntwk$to <- gsub("_HUMAN","",gsub(".*:","",ntwk$to))
head(ntwk)
seeds = "tr|Q5STU3|Q5STU3"
genes01 <- downStreamGenes(ntwk, seeds, N = 1, directed = TRUE)
genes02 <- downStreamGenes(ntwk, seeds, N = 2, directed = TRUE)
genes03 <- downStreamGenes(ntwk, seeds, N = 3, directed = TRUE)
intersect(genes01, genes03); intersect(genes01, genes02);intersect(genes02, genes03)

sig <- cbind(genes03)
net <- ntwk
net01 <- merge(net, sig, by.x = 1, by.y = 1, all = FALSE)
net02 <- merge(net01, sig, by.x = 2, by.y = 1, all = FALSE)
table(sort(union(net02$to,net02$from)) == sort(genes03))

setwd("C:/Users/wange13/Desktop/MSBB_new_rel/PHG_proteomics")
meta <- read.table("MSBB.PHG.proteomics_protein_meta.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
meta$Protein <- gsub("_HUMAN","",meta$Protein)
meta <- meta[, c(2,1,5)]
head(meta)
nodes <- meta[meta$Protein%in%genes03,]

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/bon_new/Protein_OMS")
meth <- read.table("AD_meth_signatures_proteomics.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
head(meth)
nodes$oms.neg <- nodes$Protein%in%meth[meth$module == "oms.neg",]$protein.id
nodes$oms.neg <- ifelse(nodes$oms.neg == TRUE, "negative", "none")
nodes$oms.pos <- nodes$Protein%in%meth[meth$module == "oms.pos",]$protein.id
nodes$oms.pos <- ifelse(nodes$oms.pos == TRUE, "positive", "none")
nodes$pos.neg.none <- paste(nodes$oms.neg, nodes$oms.pos, sep="_")
table(nodes$pos.neg.none)
head(nodes)

splt <- split(nodes, nodes$Symbol)
df <- NULL

for (nnn in 1:length(splt)){
  if(dim(splt[[nnn]])[1] > 1){splt[[nnn]]$Symbol <- paste(splt[[nnn]]$Symbol, (1:dim(splt[[nnn]])[1]), sep="."); df <- rbind(df, splt[[nnn]]) }
else {df <- rbind(df, splt[[nnn]])}
}
df00 <- df[,1:2]; head(df00)

nodes00 <- merge(nodes, df00, by.x = 1, by.y = 1, all.x = TRUE)
head(nodes00)
colnames(nodes00)[dim(nodes00)[2]] = "Symbol.with.dup"

setwd(wd)
dir = "bn_subnet"
if(dir.exists(dir)) {print ("This dir has already existed!!!")} else{dir.create(dir)}
seeds <- gsub("\\|","-",seeds)
write.table(net02, file= paste0(paste0("bn_subnet/PHG_proteomics_KDA_subBN_",seeds),"updated.txt"), 
            sep= "\t", row.names = FALSE, quote = FALSE)
write.table(nodes00, file= paste0(paste0("bn_subnet/PHG_proteomics_subBN_nodes_",seeds),"updated.txt"),
            sep= "\t", row.names = FALSE, quote = FALSE)




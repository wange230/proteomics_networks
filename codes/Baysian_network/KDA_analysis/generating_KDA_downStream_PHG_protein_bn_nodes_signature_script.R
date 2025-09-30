### bn networks finding downstream targets
rm(list=ls())
wd = "C:/Users/wange13/Desktop/MSBB_new_rel"
setwd(wd)
kda <- read.table("KeyDrivers/AD_L0_KDx_combined.xls", sep="\t", stringsAsFactors = FALSE, header = TRUE)
kda <- kda[,c(2,1)]
kda$keydrivers <- gsub(".*:","",kda$keydrivers)
head(kda)

df = NULL
for (nn in 1:dim(kda)[1]) {
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/bn_network")
library(KDA)
ntwk <- read.table("global.BN_cys.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)
head(ntwk)
ntwk$Source <- gsub(".*:","", ntwk$Source)
ntwk$Target <- gsub(".*:","", ntwk$Target)
colnames(ntwk) <- c("from","to")
head(ntwk)
#seeds = "sp|P26038|MOES_HUMAN"
seeds = kda$keydrivers[nn]
genes03 <- downStreamGenes(ntwk, seeds, N = 3, directed = TRUE)
net <- data.frame(protein.id = genes03, module = rep(paste0(seeds, ".Layers3"), length(genes03)))
df <- rbind(df, net)
}
head(df)

setwd(wd)
dir = "bn_subnet"
if(dir.exists(dir)) {print ("This dir has already existed!!!")} else{dir.create(dir)}
write.table(df, file= "bn_subnet/PHG_protein_KDA_BN_nodes_signatue.txt", sep= "\t", row.names = FALSE, quote = FALSE)



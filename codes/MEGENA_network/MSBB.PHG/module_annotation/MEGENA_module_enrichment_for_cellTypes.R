##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new/hyperGeTest_new/PHG_protein_mol_ranking"
setwd(wd)
library(purrr)
library(dplyr)

bigtbl <- read.table("PHG_proteomics_DEP_FC11_based_MEGENA_module_ranking_07_21_2020_withNewCellType_GO.txt",
                      sep="\t", stringsAsFactors = F, header = T)
row.names(bigtbl) = bigtbl$Module.id
df = bigtbl[c("moduleRankingOrder",
              "ast",
              "end",
              "mic",                         
              "neu",                         
              "oli" ,                        
              "opc" )]
head(df)
df1 = df[,-1]
mat = as.matrix(df1)
#apply(mat,1, which.min)
ind = apply(mat,1, which.min)
p <- c()
for(ii in 1:nrow(mat)){
p = c(p, mat[ii,][ind[ii]])
}
cts = data.frame(celltype = names(p), p.min = p)
cts$mod = row.names(mat)
head(cts)

tbl = merge(df, cts, by.x=0, by.y = 3, all = F)
head(tbl)

tbl00 = tbl[tbl$p.min<0.05,]
colnames(tbl00)[1] = "modules"
setwd(wd)
write.table(tbl00, file= "MSBB_PHG_MEGENA_module_enrichment_for_celltypes_FDR05.txt",
            sep="\t", row.names = F, quote = F)
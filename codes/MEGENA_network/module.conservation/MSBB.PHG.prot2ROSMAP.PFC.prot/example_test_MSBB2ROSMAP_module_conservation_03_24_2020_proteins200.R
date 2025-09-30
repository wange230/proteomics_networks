rm(list = ls())

simple.convert = function(x) 
{
  res = mapply(FUN = function(x,y) data.frame(id = x,module = rep(y,length(x)),stringsAsFactors = FALSE),x = x,y = as.list(names(x)),SIMPLIFY = FALSE)
  res = do.call('rbind.data.frame',res)
  rownames(res) = NULL
  return(res)
}

read.geneSet <- function(geneSet.file)
{
 gene.list <- readLines(geneSet.file)
 gene.list <- strsplit(gene.list,"\t")
 
 names(gene.list) <- sapply(gene.list,function(x) x[1])
 gene.list <- lapply(gene.list,function(x) x[2:length(x)])
 return(gene.list)
}

library(WGCNA)

####################
wd <- "/sc/hydra/projects/zhangb03a/shared/DET_WINA_MEGENA/PHG_ROSMAP_proteinMEGENA_module_conservation/MSBB_ROSMAP_conservation"
setwd(wd)
out.dir <- "module_preservation_by_WGCNA";
outfname <- "TCGA_GCC.vs.TCGA_vGCCrna9209perm200"
dir.create(out.dir)
load("RSOMAP_MSBB_proteomics_conservation.RData")
data.mat <- msbb01

valid.mat <- ros01

common.gene = intersect(rownames(data.mat),rownames(valid.mat))
multiExpr = list(TCGA = list(data = t(data.mat[match(common.gene,rownames(data.mat)),])),
                 vTCGA = list(data = t(valid.mat[match(common.gene,rownames(valid.mat)),])))

# subset modules for intersection
mol <- read.table("PHG_Proteomics_MEGENA_reformat00.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
mol00 <- mol[, c(5,4)]
colnames(mol00) <- c("id","module")
mol00$module <- gsub("c","M",mol00$module)
head(mol00)
mdf.o <- mol00
common.mdf = subset(mdf.o,id %in% common.gene)

# make duplicate labels to handle multiplicity
common.mdf$new.id = paste(common.mdf$id,"|",common.mdf$module,sep = "")
multiExpr = lapply(multiExpr,function(x,y) {out = x$data[,match(y$id,colnames(x$data))];colnames(out) = y$new.id;list(data = out)},y = common.mdf)
multiColor = list(TCGA = common.mdf$module)

max.module.size = 2000

mp = modulePreservation(multiExpr, multiColor,
                        dataIsExpr = TRUE,
                        referenceNetworks = 1,
                        nPermutations = 200,
                        randomSeed = 5,
                        maxGoldModuleSize = max.module.size+1, 
                        maxModuleSize = max.module.size + 1,
                        quickCor = 1, verbose = 3)


saveRDS(mp,file = paste(out.dir,"/",outfname,".RDS",sep = ""))

quit(save = "no",status = 0)			

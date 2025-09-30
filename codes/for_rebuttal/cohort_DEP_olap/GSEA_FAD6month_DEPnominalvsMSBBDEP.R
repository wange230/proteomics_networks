## mouse 6 month DEP vs MSBB DEP
rm(list=ls())
wd = "C:/Users/wange13/Desktop/DEP.validation/"

# prepare genelist
setwd("C:/Users/wange13/Desktop/mouse5fad/mouse_FAD_vs_WT_DEP")
mo.dep <- read.table("mouse_FAD_vs_WT_month6_Sex_adj_DEP.txt",
                     sep="\t", stringsAsFactors = F, header = T)
head(mo.dep)

# mo6dep <- mo.dep[mo.dep$adj.P.Val < 0.05 & abs(mo.dep$logFC)>log2(1.17),]
# mo6dep <- mo.dep[mo.dep$adj.P.Val < 0.05,]
mo6dep <- mo.dep[mo.dep$P.Value < 0.05,]
head(mo6dep)
mo6dep = mo6dep[c("HGNC.symbol","logFC","P.Value")]
mo6dep$module = rep("combined",nrow(mo6dep) )
mo6dep = mo6dep[!is.na(mo6dep$HGNC.symbol),]
mo6dep = mo6dep[order(mo6dep$logFC, decreasing = T),]
head(mo6dep)
genelist = mo6dep$logFC
names(genelist) = mo6dep$HGNC.symbol

# prepare term2gene
setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP")
all <- read.table("MSBB_ROSMAP_protein_signatures_at_protein_GN.level.txt",
                   sep="\t", stringsAsFactors = F, header = T)
phg <- all[all$cohort == "MSBB",]
phg = phg[!is.na(phg$protein.id),]
term2gene = phg[c("cohort","protein.id")]
head(term2gene)

#intersect(names(genelist),term2gene$protein.id)
library(clusterProfiler)
results <- GSEA(
  geneList = genelist,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 5000,
  eps = 1e-50,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE = term2gene,
  #TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea"
)
head(results)
df = results@result

# output
setwd(wd)
dir.create("gsea.results")
library(enrichplot)

pdf(file="gsea.results/FAD6mon_GSEA_for_MSBB_DEP_without_pvaluey_table.pdf",
    width = 5, height = 4.5)
gseaplot2(results,geneSetID = c(1), pvalue_table = FALSE)
dev.off()

pdf(file="gsea.results/FAD6mon_GSEA_for_MSBB_DEP_with_pvaluey_table.pdf",
     width = 5, height = 4.5)
gseaplot2(results,geneSetID = c(1), pvalue_table = TRUE)
dev.off()

write.table(df, file="gsea.results/5xFAD6mon_GSEA_for_MSBB_DEP.txt",
            sep="\t", row.names = F, quote = F)
#
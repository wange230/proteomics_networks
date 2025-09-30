##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/revision"
setwd(wd)
dir.create("Figs2024")
library(EnhancedVolcano)
library(airway)
library(magrittr)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019")
meta = read.table("protein.meta.with.uniqueSymbol.txt", sep="\t", stringsAsFactors = F, header = T)
head(meta)

cerad <- read.table("Protein_exprCorrected_HvsL_CERAD_final.tsv",
                     sep="\t", stringsAsFactors = F, header = T)
head(cerad)
cerad_meta <- merge(meta,cerad, by.x = 0, by.y =0, all =FALSE)
cerad_meta <- cerad_meta[order(cerad_meta$adj.P.Val, decreasing = F),]
row.names(cerad_meta) <- cerad_meta$protein.gene
head(cerad_meta)
summary(cerad_meta$logFC)
sel = cerad_meta[abs(cerad_meta$logFC)>log2(1.1),]
sel.symbol <- rownames(sel)[1:10]
sel.symbol = c(sel.symbol,"TTBK2.2","S100A10.1","NPTX2.1","NRN1.1","VGF.1","RPH3A.2","OLFM3.1","SYT12.1
")


table(row.names(cerad_meta) == rownames(cerad_meta))
summary(-log10(cerad_meta$adj.P.Val))
p5 <- EnhancedVolcano(cerad_meta,
                      lab = rownames(cerad_meta),
                      x = "logFC",
                      y = "adj.P.Val",
                      drawConnectors = TRUE,
                      lengthConnectors = unit(0.005, "npc"),
                      selectLab = sel.symbol, 
                      xlab = bquote(~Log[2]~ "(FC)"),
                      ylab = bquote(~-Log[10]~italic((p.adj))),
                      pCutoff = 0.05,
                      FCcutoff = log2(1.1),
                      cutoffLineWidth = 0.5,
                      xlim = c(summary(cerad_meta$logFC)[1]-0.25,
                               summary(cerad_meta$logFC)[6]+0.25),
                      ylim = c(0.0, summary(-log10(cerad_meta$adj.P.Val))[6]+2),
                      labSize = 3.0,
                      pointSize = 1.0,
                      axisLabSize = 10,
                      colAlpha = 0.5,
                      legendLabels = c("NS",
                                       "log2(FC)",
                                       "p.adj",
                                       "p.adj & log2(FC)"),
                      legendPosition = "bottom",
                      legendLabSize = 10,
                      legendIconSize = 2.5,
                      title = " ", 
                      subtitle = "DefiniteAD vs. NL",
                      subtitleLabSize = 12.5,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      borderWidth = 0.75
                      )
#p5
library(gridExtra)
library(grid)
setwd(wd)
pdf(file = "Figs2024/DEP.PHG.CERAD00.pdf", 
      width = 5.5, height = 5.5)
print(p5)
dev.off()

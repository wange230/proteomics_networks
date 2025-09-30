##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/ahnak_proteomics/DEP"
setwd(wd)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(gridExtra)
library(grid)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
organism = org.Mm.eg.db
## plot 
df = read.table("AHNAK_KD_vs_CTRL_normalization.tsv",sep="\t", 
                  stringsAsFactors = F, header = T)
head(df)
splt = split(df,df$GN)
dff <- NULL
for(nn in 1:length(splt)){
dat = splt[[nn]]
dat = dat[order(dat$P.Value),]
dat = dat[1,]
dff = rbind(dff,dat)
}
df = dff
row.names(df) = df$GN
df = df[order(df$P.Value),]
head(df)
sel00 = c(row.names(df[df$logFC<0,])[1:10],"AHNAK")
sel01 = c(row.names(df[df$logFC>0,])[1:10])
sel = c(sel00,sel01)
sel = sel[!sel==""]
## plot 

peak <- df
head(peak)
summary(peak$logFC)
sel.symbol <- sel
summary(-log10(peak$adj.P.Val))
summary(-log10(peak$P.Value))
cutoff = 0.05

p5 <- EnhancedVolcano(peak,
                      lab = rownames(peak),
                      x = "logFC",
                      y = "P.Value",
                      drawConnectors = TRUE,
                      selectLab = sel.symbol, 
                      xlab = bquote(~Log[2]~ (FC)),
                      ylab = bquote(~-Log[10]~italic((P.Value))),
                      pCutoff = cutoff, # significant p upon adjustment
                      FCcutoff = log2(1.1),
                      cutoffLineWidth = 0.35,
                      lengthConnectors = unit(0.0, "cm"),
                      widthConnectors = 0.35,
                      xlim = c(-0.25+summary(peak$logFC)[1],summary(peak$logFC)[6]+0.25),
                      ylim = c(0.0-1, summary(-log10(peak$adj.P.Val))[6]+3),
                      labSize = 3.0,
                      pointSize = 1.0,
                      axisLabSize = 10,
                      colAlpha = 0.5,
                      legendLabels = c("NS","Log2 FC","p_value","p-value & Log2 FC"),
                      legendPosition = "bottom",
                      legendLabSize = 10,
                      legendIconSize = 5.0,
                      title = " ", subtitle = "AHNAK iPSC DEP analysis",
                      subtitleLabSize = 10.0,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      borderWidth = 0.5)
#p5
setwd(wd)
pdf (file = "AHNAK_iPSC_DEP_12_19_2024.pdf", width = 5.0, height = 6.0 )
print(p5)
dev.off()
# over
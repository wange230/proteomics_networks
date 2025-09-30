#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/signal_maps/updated"
setwd(wd)
library(dplyr)
library(purrr)
library(ggplot2)


df = read.table("AHNAK_DEP03_for_GSEA_BP_removing_go_no_annot_11012023.txt",
                 sep="\t", stringsAsFactors = F,header = T)
head(df)
df$`-log10(qvalue)` = log10(df$qvalue.child) * (-1)

summary(df$`-log10(qvalue)`)
g = ggplot(data = df, aes(x = NES.child, y = qvalue.child, size = `-log10(qvalue)`)) +
    geom_point(shape=15, color="gray10")
g0 = g + theme(legend.position="right",
               legend.title = element_text(color = "blue", size = 20, face ="bold"),
               legend.text = element_text(color = "blue", size = 15))

pdf(file="legend.pdf")
print(g0)
dev.off()
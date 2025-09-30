##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP"
dir.create("DEP_DEG_diff")
setwd(wd)
setwd("C:/Users/wange13/Desktop/ROSMAP_meth/DEP")
deg <- read.table("ROSMAP_consensus_DEG.txt",
                  sep="\t", stringsAsFactors = F, header = T)
head(deg)
length(deg$Symbol[deg$comp_DEG == "Up"])+length(deg$Symbol[deg$comp_DEG == "Dn"])
setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP")
dep <- read.table("ROSMAP_consensus_DEP.txt",
                  sep="\t", stringsAsFactors = F, header = T)
head(dep)
table(dep$symbol%in%deg$Symbol)

#superexactTest
background.signature = intersect(deg$Symbol,dep$symbol)
length(background.signature)

dep.sel <- dep[dep$symbol%in%background.signature,]

deg.sel <- deg[deg$Symbol%in%background.signature,]

#superExactTest
a=as.character(dep.sel$symbol[dep.sel$comp_DEP == "Up"]) 
b=as.character(dep.sel$symbol[dep.sel$comp_DEP == "Dn"])
c=as.character(dep.sel$symbol[dep.sel$comp_DEP == "Neutral"])
d=as.character(deg.sel$Symbol[deg.sel$comp_DEG == "Up"]) 
e=as.character(deg.sel$Symbol[deg.sel$comp_DEG == "Dn"])
f=as.character(deg.sel$Symbol[deg.sel$comp_DEG == "Neutral"])

library(ggVennDiagram)
library(ggplot2)
library("ggsci")
x <- list(`DEP(Up)` = a, `DEG(Up)` = d)
pdf(file="ROSMAP_DEP_upvsDEG_up2024.pdf",height = 3.5, width = 3.5)
p= ggVennDiagram(x, label_color = "black",label_alpha = 0,
                 set_size = 5.5,
                 label_size = 6.0,
                 label = "count",
                 set_color = "blue",
                 label_txtWidth = 80,
                 #edge_lty = "dashed",
                 #edge_size = 1,
                 show.legend = FALSE
                 #edge_color = "white"
)
p + scale_fill_distiller(palette = "Pastel1") +
  scale_color_manual(values = c("blue","blue")) +
  theme(legend.position = "none")
dev.off()


x <- list(`DEP(Down)` = b, `DEG(Down)` = e)
pdf(file="ROSMAP_DEP_dnvsDEG_dn2024.pdf",height = 3.5, width = 3.5)
p= ggVennDiagram(x, label_color = "black",label_alpha = 0,
                 set_size = 5.5,
                 label_size = 5.0,
                 label = "count",
                 set_color = "blue",
                 label_txtWidth = 80,
                 #edge_lty = "dashed",
                 #edge_size = 1,
                 show.legend = FALSE
                 #edge_color = "white"
)
p + scale_fill_distiller(palette = "Pastel1") +
  scale_color_manual(values = c("blue","blue")) +
  theme(legend.position = "none")
dev.off()

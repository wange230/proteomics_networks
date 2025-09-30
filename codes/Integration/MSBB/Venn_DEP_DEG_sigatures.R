##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/volcano_plot"
setwd(wd)
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/MSBB_RNAseq_latest/mixed_model_correction_gene/comp_DEG")
deg <- read.table("summay_DEG_across_trait_in_BM36_region_withVotes_usingGeneSymbolUnique_updated00.txt",
                  sep="\t", stringsAsFactors = F, header = T)
head(deg)
length(deg$GenSymbol[deg$comp.ADvsNL == "Up"])+length(deg$GenSymbol[deg$comp.ADvsNL == "Dn"])
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019/DEP_union")
dep <- read.table("Protein_exprCorrected_HvsL_Braak_CDR_CERAD_PLQ_withVotes_uniqueGenSymbol_updated.txt",
                  sep="\t", stringsAsFactors = F, header = T)
head(dep)
table(dep$Symbol%in%deg$GenSymbol)
dep.only <- dep[!dep$Symbol%in%deg$GenSymbol,]
dep.only.sel <- dep.only[,c(2,dim(dep.only)[2])]
dep.only.sel$comp.ADvsNL = paste("protein_sepcific",dep.only.sel$comp.ADvsNL, sep=".")
colnames(dep.only.sel)[1] <- "GN"
head(dep.only.sel)
table(deg$GenSymbol%in%dep$Symbol)
deg.only <- deg[!deg$GenSymbol%in%dep$Symbol,]
deg.only.sel <- deg.only[,c(2,dim(deg.only)[2])]
deg.only.sel$comp.ADvsNL = paste("mRNA_sepcific",deg.only.sel$comp.ADvsNL, sep=".")
colnames(deg.only.sel)[1] <- "GN"
head(deg.only.sel)

tbl <- rbind(dep.only.sel,deg.only.sel)
head(tbl)
length(unique(tbl$GN))
setwd(wd)
dir = "DEP_DEG_diff"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
#write.table(dep.only, file = 'DEP_DEG_diff/PHG_DEP_only_signatures_updated.txt',
 #           row.names = FALSE, quote = F, sep = "\t")
#write.table(deg.only, file = 'DEP_DEG_diff/PHG_DEG_only_signatures_updated.txt', 
 #           row.names = FALSE, quote = F, sep = "\t")
#write.table(tbl, file = 'DEP_DEG_diff/PHG_protein_or_mRNA_specific_signatures_updated.txt', 
 #           row.names = FALSE, quote = F, sep = "\t")

#superexactTest
background.signature = intersect(deg$GenSymbol,dep$Symbol)
length(background.signature)

dep.sel <- dep[dep$Symbol%in%background.signature,]
intersect(dep.sel$Symbol,tbl$GN) # should 0
deg.sel <- deg[deg$GenSymbol%in%background.signature,]
intersect(deg.sel$Symbol,tbl$GN) # should 0
#superExactTest
a=as.character(dep.sel$Symbol[dep.sel$comp.ADvsNL == "Up"]) # table(dep.sel$comp.ADvsNL)
b=as.character(dep.sel$Symbol[dep.sel$comp.ADvsNL == "Down"])
c=as.character(dep.sel$Symbol[dep.sel$comp.ADvsNL == "Neutral"])
d=as.character(deg.sel$GenSymbol[deg.sel$comp.ADvsNL == "Up"]) ## table(deg.sel$comp.ADvsNL)
e=as.character(deg.sel$GenSymbol[deg.sel$comp.ADvsNL == "Dn"])
f=as.character(deg.sel$GenSymbol[deg.sel$comp.ADvsNL == "Neut"])

library(ggVennDiagram)
library(ggplot2)
library("ggsci")
x <- list(`DEP(Up)` = a, `DEG(Up)` = d)
pdf(file="MSBB_DEP_upvsDEG_up2024.pdf",height = 3.5, width = 3.5)
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
pdf(file="MSBB_DEP_dnvsDEG_dn2024.pdf",height = 3.5, width = 3.5)
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


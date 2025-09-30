##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/volcano_plot"
setwd(wd)
library(EnhancedVolcano)
library(airway)
library(magrittr)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019")
meta = read.table("protein.meta.with.uniqueSymbol.txt", sep="\t", stringsAsFactors = F, header = T)
head(meta)


cdr <- read.table("Protein_exprCorrected_HvsL_CDR_final.tsv", 
                  sep="\t", stringsAsFactors = F, header = T)
head(cdr)
cdr_meta <- merge(meta,cdr, by.x = 0, by.y =0, all =FALSE)
cdr_meta <- cdr_meta[order(cdr_meta$adj.P.Val, decreasing = F),]
row.names(cdr_meta) <- cdr_meta$protein.gene
head(cdr_meta)
dep.cdr = cdr_meta[abs(cdr_meta$logFC)>log2(1.1)&cdr_meta$adj.P.Val<0.05,]
dep.cdr$sign = ifelse(dep.cdr$logFC>0, "Up","Dn")


cerad <- read.table("Protein_exprCorrected_HvsL_CERAD_final.tsv", 
                    sep="\t", stringsAsFactors = F, header = T)
head(cerad)
cerad_meta <- merge(meta,cerad, by.x = 0, by.y =0, all =FALSE)
cerad_meta <- cerad_meta[order(cerad_meta$adj.P.Val, decreasing = F),]
row.names(cerad_meta) <- cerad_meta$protein.gene
head(cerad_meta)
dep.cerad = cerad_meta[abs(cerad_meta$logFC)>log2(1.1)&cerad_meta$adj.P.Val<0.05,]
dep.cerad$sign = ifelse(dep.cerad$logFC>0, "Up","Dn")

a1 = dep.cdr$Protein[dep.cdr$sign=="Dn"]
a2 = dep.cdr$Protein[dep.cdr$sign=="Up"]
aa1 = dep.cerad$Protein[dep.cerad$sign=="Dn"]
aa2 = dep.cerad$Protein[dep.cerad$sign=="Up"]

source("C:/pc.backup/Documents/R_code/R_functions_for_cor_network/Function_4signature_venn_pdf.R")
#fix(venn4sig)
setwd(wd)
dir.create("CDRvsCEARD")
setwd("CDRvsCEARD")

#venn4sig(a1,a2,aa1,aa2)

# FET test
bkg = cdr_meta$Protein ## Background size
bkg.gen = unique(cdr_meta$Symbol[!is.na(cdr_meta$Symbol)])
# positive set
dat = dep.cdr
dat = dat[dat$Protein%in%bkg,]
head(dat)
set1 <- dat[,c(3,ncol(dat))]
head(set1)
splt1 <- split(set1, set1$sign)
length(splt1)

## sampling set
datt= dep.cerad
datt = datt[datt$Protein%in%bkg,]
head(datt)
set2 <- datt[,c(3,ncol(datt))]
head(set2)
splt2 <- split(set2, set2$sign)
length(splt2)

source("c:/pc.backup/Documents/R_code/R_functions_for_cor_network/geneSetOverLap.R")
size = length(bkg)  # 12147

hyperGeTest <- GeneSetOverlap(set1,set2, Background = "large",BackgroundSize = size)
enTable <- hyperGeTest$EnrichTable
head(enTable)
enTable$p.value.FDR <- p.adjust(enTable$`P value`, method = "fdr", n = length(splt1)*length(splt2) )
enTable <- enTable[, c(1:8, 10,9)]

setwd(wd)
setwd("CDRvsCEARD")
dir ="hyperGeTest"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
#write.table(enTable, file = "hyperGeTest/EnrichmentTable_DEP_CDRvsDEP_CERAD_MSBB.txt",
 #          sep="\t", row.names = FALSE)

# GO anaysis
setwd(wd)
setwd("CDRvsCEARD")
dir.create("GO")
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(clusterProfiler)
organism = org.Hs.eg.db

CDR_sep_dn = unique(cdr_meta$Symbol[cdr_meta$Protein%in%setdiff(a1, aa1)])
CDR_sep_up = unique(cdr_meta$Symbol[cdr_meta$Protein%in%setdiff(a2, aa2)])
CDR_CERAD_dn = unique(cdr_meta$Symbol[cdr_meta$Protein%in%intersect(a1, aa1)])
CDR_CERAD_up = unique(cdr_meta$Symbol[cdr_meta$Protein%in%intersect(a2, aa2)])
CERAD_sep_dn = unique(cdr_meta$Symbol[cdr_meta$Protein%in%setdiff(aa1, a1)])
CERAD_sep_up = unique(cdr_meta$Symbol[cdr_meta$Protein%in%setdiff(aa2, a2)])
lt = list(CDR_sep_dn=CDR_sep_dn,
          CDR_sep_up=CDR_sep_up,
          CDR_olap_CERAD_dn=CDR_CERAD_dn,
          CDR_olap_CERAD_up=CDR_CERAD_up,
          CERAD_sep_dn=CERAD_sep_dn,
          CERAD_sep_up=CERAD_sep_up)
nams = names(lt)
for(ii in 1:length(nams)){
#ii =1
gene = lt[[nams[ii]]]
ego <- enrichGO(gene          = gene,
                keyType = "SYMBOL",
                universe      = bkg.gen,
                OrgDb         = organism,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                minGSSize = 10, 
                maxGSSize = 5000,
                readable      = TRUE)
head(ego)
results = ego@result  

if(nrow(results)<1){next}
write.table(results, file= paste0("GO/enrichGO_",nams[ii],"_9209bkg.txt"),
            sep="\t", row.names = F, quote = F)
rm(gene, ego, results)
}# ii


# legend
rm(list=ls())
wd = "C:/Users/wange13/Desktop/volcano_plot/CDRvsCEARD"
setwd(wd)
#
library(ComplexHeatmap)

lgd = Legend(labels = c("CDR.down","CDR.up","CERAD.down","CERAD.up"), title = "", 
             legend_gp = gpar(fill = c("orange", "red", "turquoise", "blue")),
             labels_gp = gpar(col = "black", fontsize = 8)
)

pdf("legends.pdf", width = 1.0, height = 1.0)
draw(lgd)
dev.off()


# basic R
barplot(table(mtcars$gear), col = 1:4)

legend("topright",
       legend = c("CDR.down","CDR.up","CERAD.down","CERAD.up"),
       fill = c("orange", "red", "turquoise", "blue"),       # Color of the squares
       border = "white")
       #border = "black") # Color of the border of the squares


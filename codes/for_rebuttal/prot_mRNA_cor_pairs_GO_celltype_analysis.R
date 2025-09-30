##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/revision/revision_R1/GO2024"

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019/Protein_mRNA_cor_forDGCA")
ccc = read.table("PHG_protein__cor_mRNA_uniqueSymbol.8498.txt",
                  sep="\t", stringsAsFactors = F, header = T)

head(ccc)
table(is.na(ccc$Gene.Symbol))

sig = ccc[ccc$p.adjust.BH<0.05,]
table(sig$r.value>0)
summary(abs(sig$r.value))
sig = sig[order(abs(sig$r.value), decreasing = T),]
head(sig)
sig.r051 = sig[abs(sig$r.value)>0.5,]
table(sig.r051$r.value>0)

sig.r050 = sig[!abs(sig$r.value)>0.5,]

insig = ccc[!ccc$p.adjust.BH<0.05,]
insig = insig[order(abs(insig$r.value), decreasing = T),]
head(insig)
table(insig$r.value>0)
insig.bottom200 = insig[(nrow(insig)-199):nrow(insig),]
tail(insig.bottom200)

setwd("C:/Users/wange13/Desktop/volcano_plot/DEP_DEG_diff/DEP_DEG.signature")

dep.deg = read.table("PHG_DEG_DEP_combination_signatures.tsv", sep="\t", stringsAsFactors = F, header = T)

dep.deg.up = dep.deg$GN[dep.deg$module=="DEP(up) & DEG(up)"]
dep.deg.dn = dep.deg$GN[dep.deg$module=="DEP(down) & DEG(down)"]

length(intersect(dep.deg.up, sig$Gene.Symbol)) # 503/563 #0.8934281
length(intersect(dep.deg.dn, sig$Gene.Symbol)) # 265/278 #0.9532374 
# together (503+265)/(563+278) = 0.9131986

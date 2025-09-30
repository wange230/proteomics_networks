##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP.updated"
setwd(wd)
library(dplyr)
library(purrr)

fls = list.files()
fls = fls[grepl("\\.tsv", fls)]
head(fls)
contrasts = gsub(".*Corrected_","", fls)
contrasts = gsub(".*corrected_","", contrasts)
traits = gsub("\\.tsv","",gsub(".*_","",contrasts))
contrasts = gsub("_.*","", contrasts)
contrasts = gsub("HvsL","High.AD.pathology_vs_Low.AD.pathology",contrasts)
contrasts = gsub("HvsM","High.AD.pathology_vs_Moderate.AD.pathology",contrasts)
contrasts = gsub("MvsL","Moderate.AD.pathology_vs_Low.AD.pathology",contrasts)

#fls
#traits
#contrasts

df <- NULL
for(ii in 1:length(fls)){
# ii =1
dat = read.table(fls[ii], sep = "\t", stringsAsFactors = F, header = T)
dat$ChangeDirection = ifelse(dat$logFC>0, "Up", "Down")
dat$Contrast = rep(contrasts[ii], nrow(dat))
dat$Trait = rep(traits[ii], nrow(dat))
dat$Cohort = rep("ROSMAP", nrow(dat))
head(dat)
dat05 = dat[dat$adj.P.Val<0.05,]
if(nrow(dat)<1){next}
df = rbind(df,dat05)
rm(dat)
}
head(df)

splt = strsplit(row.names(df), "\\|")
# splt[[1]]
GN = map_chr(splt, function(x){x[[1]]})
AccessionID = map_chr(splt, function(x){x[[2]]})

df$ProteinID = row.names(df)
df$Gene.Symbol = GN
df$Protein.accession = AccessionID
#colnames(df)

dff = df[c("ProteinID",
           "Gene.Symbol",
           "Protein.accession",
           "logFC",
           "ChangeDirection",
           "P.Value",
           "adj.P.Val",
           "Contrast",
           "Trait",
           "Cohort" )]
head(dff)
setwd(wd)

write.table(dff, file="ROSMAP_DEP_summary_updated.txt",sep="\t", row.names = F, quote = T)
# over
## MSBB prot over MSBB mRNA and ROSMAP mRNA and prot MEGENA module overlap
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/output/imputed/MEGENA/FET.mod_olap"
setwd(wd)
library(tidyverse)

fls = list.files()
fl  = fls[!grepl("reformat",fls)] 
fn  = gsub("\\.txt","", fl)

for(ii in 1:length(fn)){
target = read.table(fl[ii],sep="\t",stringsAsFactors = F, header = T)
head(target)
df = target[,c(1,2,9)]
head(df)

dff = spread(df, key = Set2, value = p.value.FDR)
row.names(dff) = dff$Set1
dff = dff[,-1]
dff[dff > 0.05] = ""
dff$src.module = row.names(dff)
dff = dff[,c(ncol(dff), 1:(ncol(dff)-1))]
write.table(dff, file= paste0(fn[ii], "_reformat.txt"),
            sep="\t", row.names = F, quote = F)
rm(dff)
}# ii
##
rm(list=ls())
wd = "/sc/arion/projects/adineto/permutation_2024/"
setwd(wd)

library(purrr)
library(tidyverse)

fls = list.files("DEP_downsampling")
fn  = fls[grepl("final_",fls)]
df <- NULL
for (ii in 1:length(fn)){
# ii = 1
dat = read.table(paste0("DEP_downsampling/",fn[ii]), 
                 sep="\t", stringsAsFactors = F, header = T)
datt = dat[abs(dat$logFC)>log2(1.1)&dat$adj.P.Val<0.05,] 
datt$perm = rep(paste0("permutation_",ii), nrow(datt)  )
head(datt)
datt = datt[c("Protein","Symbol","logFC","adj.P.Val","perm")]
df = rbind(df, datt)
rm(dat, datt)
print(ii)
} # ii
head(df)
#write.table(df, file="Protein_exprCorrected_defAD_vs_NL_Female_CERAD_final_perm1000",
 #           sep="\t", row.names = F, quote = F)

df = read.table("Protein_exprCorrected_defAD_vs_NL_Female_CERAD_final_perm1000", 
                sep="\t",stringsAsFactors = F, header = T)
splt = split(df, df$Protein)
len = map_int(splt, function(x){nrow(x)})

tbl = data.frame(protein= names(len), occurrance = len)
head(tbl)
tbl200 = tbl[tbl$occurrance>200,]
tbl300 = tbl[tbl$occurrance>300,]
tbl400 = tbl[tbl$occurrance>400,]
tbl500 = tbl[tbl$occurrance>500,]
tbl600 = tbl[tbl$occurrance>600,]
tbl700 = tbl[tbl$occurrance>700,]
tbl800 = tbl[tbl$occurrance>800,]
tbl900 = tbl[tbl$occurrance>900,]

dat = data.frame(cutoff = c(200, 300,400,500,600,700,800,900),
                 DEPs = c(nrow(tbl200)
                          ,nrow(tbl300)
                          ,nrow(tbl400)
                          ,nrow(tbl500)
                          ,nrow(tbl600)
                          ,nrow(tbl700)
                          ,nrow(tbl800),nrow(tbl900)))
head(dat)
library(ggplot2)
dat$cutoff = factor(dat$cutoff)
g = ggplot(data=dat, aes(cutoff, DEPs, color = cutoff) ) +
    geom_point()
g

pdf(file="summry_cutoff_downsample_Female.pdf")
g
dev.off()
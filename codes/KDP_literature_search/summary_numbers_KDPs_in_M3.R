# RCG_literature_all stored literature search results for KDPs in module M3
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/RCG_literature_all"
setwd(wd)

fls = list.files()

lite = c()
nams = c() 
for (ii in 1:length(fls)){
#  ii = 1 
  
  if(length(readLines(fls[ii])) > 1 ){
  dat = read.table(fls[ii], sep="\t", stringsAsFactors = F, header = T,quote = '"', fill = T )
  len = length(unique(dat$PMID))
  nams = c(nams, fls[ii])
  lite = c(lite, len)
  } else {
    nams = c(nams, fls[ii])
    lite = c(lite, 0)
  }
}
df = data.frame(fl = nams, length = lite )
df$fn = gsub("_literature.*", "", gsub("Resilience_RCG_AD_","", df$fl))
head(df)
df$length

hist(df$length, breaks = 500)
summary(df$length)
table(df$length>0)
#FALSE  TRUE 
#21    87 
table(df$length>10)
#FALSE  TRUE 
# 82    26 
table(df$length>20)
#FALSE  TRUE 
#93    15 
table(df$length>30)
#FALSE  TRUE 
#98    10 

table(df$length>40)
#FALSE  TRUE 
# 99     9 
table(df$length>50)
#FALSE  TRUE 
#100     8 
table(df$length>100)
#FALSE  TRUE 
#100     8

setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/Supplementary Data/excel.files")

library(readxl)
dataset <- read_excel("supplementary.data 5.xlsx", sheet = 2, skip = 4)
dataset = data.frame(dataset)
head(dataset)
dat = dataset[!dataset$Gene.Symbol == "NA",c("Gene.Symbol","KDP.ranking.Order")]

splt = split(dat, dat$Gene.Symbol)

tbl <- NULL
for(jj in 1:length(splt)){
# jj = 1
  datt = splt[[jj]]
  if(nrow(datt) >1) {
  datt = datt[order(datt$KDP.ranking.Order, decreasing = F),]
  tbl = rbind(tbl,datt[1,])
  
  } else {
    tbl = rbind(tbl, datt)
  }
  rm(datt)
}
dim(tbl)
head(tbl)
head(df)

dff = merge(df, tbl, by.x = "fn", by.y = "Gene.Symbol", all.x = F)
head(dff)
dff = dff[order(dff$KDP.ranking.Order, decreasing = F),]
setwd(wd)
#write.table(dff, file = "../KDP_in_M3_ranking_literature.txt",
 #           sep="\t", row.names = F, quote = F)

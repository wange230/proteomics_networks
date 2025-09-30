#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP"
setwd(wd)
library(reshape)
library(readxl)
library(tidyverse)

dep = read.table("Excluded_NA50_RSOMAP_PFC_proteomics_sex.age_corrected_HvsL_Braak.tsv",
                 sep ="\t",stringsAsFactors = F, header = T)
prot = row.names(dep)
head(prot)
lt = strsplit(prot,"\\|")
symbol = unlist(map(lt, function(x)x[[1]][1]))
accession = unlist(map(lt, function(x)x[[2]][1]))
gen = data.frame(proteinID = prot, 
                 symbol=symbol,
                 accession = accession)
head(gen)

contr = c("HvsL_Braak","HvsL_CERAD","HvsL_MMSE")
trait = c("Braak","CERAD","MMSE")
df = NULL
for (ii in 1:length(contr)){
dat = read.table(paste0("Excluded_NA50_RSOMAP_PFC_proteomics_sex.age_Corrected_",contr[ii],".tsv"),
                  sep="\t",stringsAsFactors = F, header = T)
head(dat)
dat$Geneid = row.names(dat)

head(dat)
dat00 = dat[c("Geneid","logFC","adj.P.Val")]
gen00 = gen
gen00 = merge(gen00, dat00, by.x = 1, by.y = 1, all = F)
head(gen00)
gen00 = gen00[!gen00$symbol == 0,]
gen00 = gen00[!is.na(gen00$symbol),]

splt = split(gen00, gen00$symbol)
dattt <-NULL
for(jj in 1:length(splt)){
datt = splt[[jj]]
if(nrow(datt) == 1) {dattt = rbind(dattt,datt)} else {
datt = datt[order(datt$adj.P.Val),]; dattt = rbind(dattt, datt[1,])}
} # jj
dattt$DEP = ifelse(dattt$logFC>0 & dattt$adj.P.Val<0.05, "Up",
                 ifelse(dattt$logFC <0 & dattt$adj.P.Val<0.05, "Dn", "Neutral"))
dattt$trait = rep(trait[ii], nrow(dattt))
df = rbind(df, dattt)
rm(gen00, datt,dattt, splt)
} # ii
head(df)

splt2 = split(df, df$trait)
names(splt2)
rm(gen)
gen <- splt2[[1]][1:3]
for(nn in 1:length(splt2)){
gen = merge(gen, splt2[[nn]][c("proteinID","DEP")], by.x =1, by.y =1, all =F)  
colnames(gen)[ncol(gen)] = names(splt2)[nn]
}
head(gen)

table(gen$Braak)
table(gen$CERAD)
table(gen$MMSE)

gen$Braak.sign = ifelse(gen$Braak=="Neutral",0, 
                        ifelse(gen$Braak=="Up",1,-1))

gen$CERAD.sign = ifelse(gen$CERAD=="Neutral",0, 
                        ifelse(gen$CERAD=="Up",1,-1))

gen$MMSE.sign = ifelse(gen$MMSE=="Neutral",0, 
                        ifelse(gen$MMSE=="Up",1,-1))

table(gen$Braak.sign)
table(gen$CERAD.sign)
table(gen$MMSE.sign)

row.names(gen) = gen$proteinID
cts = gen[c("Braak.sign","CERAD.sign","MMSE.sign")]
head(cts)
count.sum = rowSums (cts)
table(count.sum)
count0 =  rowSums(cts == 0)  
count1 =  rowSums(cts == 1)  
count_1 = rowSums(cts == -1)  
cts$count.sum = count.sum
cts$count0 = count0
cts$count1 = count1
cts$count_1 = count_1
head(cts)

cts$comp_DEP = ifelse(count1 == 3 | count1 == 2 | count1 == 1 & count_1 == 0,"Up",
                  ifelse(count_1 == 3 | count_1 == 2 | count_1 == 1 & count1 == 0,"Dn","Neutral"))
table(cts$comp_DEP)
cts00 = cts[c("count0", "count1","count_1","count.sum","comp_DEP")]
tbl = merge(gen, cts00, by.x = 1, by.y = 0, all = F)
head(tbl)
table(tbl$comp_DEP)

setwd(wd)
write.table(tbl, file="ROSMAP_consensus_DEP.txt", sep="\t", row.names = F, quote = F )
# over

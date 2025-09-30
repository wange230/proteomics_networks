##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019"
setwd(wd)

df <- read.table("DEP_union/Protein_exprCorrected_HvsL_Braak_CDR_CERAD_PLQ_final.tsv",
                  sep="\t", stringsAsFactors = F, header = T )
head(df)
row.names(df) <- df$Protein.id
df$Braak <- rep(0, dim(df)[1])
df$CDR <- rep(0, dim(df)[1])
df$CERAD <- rep(0, dim(df)[1])
df$PLQ <- rep(0, dim(df)[1])
head(df)

for(nn in 1:dim(df)[1]) {
if(df$Braak.logFC[nn] > log2(1.1) & df$Braak.adj.P.Val[nn] <0.05) {df$Braak[nn] = 1} 
if(df$Braak.logFC[nn] < -log2(1.1) & df$Braak.adj.P.Val[nn] <0.05) {df$Braak[nn] = -1}
}
head(df)

for(nn in 1:dim(df)[1]) {
  if(df$CDR.logFC[nn] > log2(1.1) & df$CDR.adj.P.Val[nn] <0.05) {df$CDR[nn] = 1} 
  if(df$CDR.logFC[nn] < -log2(1.1) & df$CDR.adj.P.Val[nn] <0.05) {df$CDR[nn] = -1}
}
head(df)

for(nn in 1:dim(df)[1]) {
  if(df$CERAD.logFC[nn] > log2(1.1) & df$CERAD.adj.P.Val[nn] <0.05) {df$CERAD[nn] = 1} 
  if(df$CERAD.logFC[nn] < -log2(1.1) & df$CERAD.adj.P.Val[nn] <0.05) {df$CERAD[nn] = -1}
}
head(df)

for(nn in 1:dim(df)[1]) {
  if(df$PLQ.logFC[nn] > log2(1.1) & df$PLQ.adj.P.Val[nn] <0.05) {df$PLQ[nn] = 1} 
  if(df$PLQ.logFC[nn] < -log2(1.1) & df$PLQ.adj.P.Val[nn] <0.05) {df$PLQ[nn] = -1}
}
head(df)

df$ranking.sum <- rep("", dim(df)[1])

for (ll in 1:dim(df)[1]){
df$ranking.sum[ll] = df$Braak[ll]+df$CDR[ll]+df$CERAD[ll]+df$PLQ[ll]
}
head(df)
table(df$ranking.sum)

##
dat <- df[,10:13]
head(dat)

keep0 <- dat == 0
sum0 <- rowSums(keep0)
df$count0 <- sum0
head(df)

keep1 <- dat == 1
sum1 <- rowSums(keep1)
df$count1 <- sum1
head(df)

keep_1 <- dat == -1
sum_1 <- rowSums(keep_1)
df$count_1 <- sum_1
head(df)
table(df$count0)
table(df$count1)
table(df$count_1)
##
df.count0 = df[df$count0 == 4,]
df.count1 = df[df$count1 > 2,]
df.count_1 = df[df$count_1 > 2,]

df.count2p = df[df$count1 == 2 & df$count_1 <2, ]
df.count2n = df[df$count1 < 2 & df$count_1 == 2, ]

df.count1p = df[df$count1 == 1 & df$count_1 < 1, ]
df.count1n = df[df$count1 < 1 & df$count_1 == 1, ]

df.pos <- rbind(rbind(df.count1,df.count2p),df.count1p)
df.nega <- rbind(rbind(df.count_1,df.count2n),df.count1n)
df.neut <- df.count0
# 9884 + 1616 + 647
df.pos$comp.ADvsNL = rep("Up", dim(df.pos)[1])
df.nega$comp.ADvsNL = rep("Down", dim(df.nega)[1])
df.neut$comp.ADvsNL = rep("Neutral", dim(df.neut)[1])

final <- rbind(rbind(df.nega,df.pos),df.neut)
head(final)
#write.table(final, file = "DEP_union/Protein_exprCorrected_HvsL_Braak_CDR_CERAD_PLQ_withVotes_new.tsv",
 #            sep="\t", row.names = F, quote = F)

## getwd()
pro.meta <- read.table("MSBB.PHG.proteomics_protein_meta_final.txt",
                        sep="\t", stringsAsFactors = F, header = T)
head(pro.meta)
final00 <- merge(final, pro.meta, by.x = 1, by.y = 1, all = F)
head(final00)
final00.na <- final00[is.na(final00$Symbol),] ## 63 with no symbol
final00 <- final00[!is.na(final00$Symbol),] ## 12084 with gene symbol

## dedupping over gene symbol across groups
table(final00$comp.ADvsNL)
nega = final00[final00$comp.ADvsNL == "Down",]
pos = final00[final00$comp.ADvsNL == "Up",]
neutr = final00[final00$comp.ADvsNL == "Neutral",]
intersect(nega$Symbol,pos$Symbol)
intersect(nega$Symbol,neutr$Symbol)
intersect(pos$Symbol,neutr$Symbol)
drop00 = intersect(nega$Symbol,pos$Symbol) ## "STXBP2"
drops01 = union(intersect(nega$Symbol,neutr$Symbol),intersect(pos$Symbol,neutr$Symbol)) ## 136

nega.sel = nega[!nega$Symbol%in%drop00,]
pos.sel = pos[!pos$Symbol%in%drop00,]
neutr.sel <- neutr[!neutr$Symbol%in%drops01,]
intersect(nega.sel$Symbol,pos.sel$Symbol) ## should 0
intersect(nega.sel$Symbol,neutr.sel$Symbol) ## should 0
intersect(pos.sel$Symbol,neutr.sel$Symbol) ## should 0
nega.drop00 = nega[nega$Symbol%in%drop00,]
pos.drop00 = pos[pos$Symbol%in%drop00,]
df.drop00 = rbind(nega.drop00,pos.drop00)
df.drop00$comp.ADvsNL <- rep("Neutral",dim(df.drop00)[1])

final01 <- rbind(rbind(rbind(nega.sel,pos.sel),neutr.sel),df.drop00)
## dedupping over gene symbol within groups
all <- NULL
splt <- split(final01, final01$comp.ADvsNL)
for (mm in 1:length(splt)) {
  ttta<- splt[[mm]]
  head(ttta); dim(ttta)
  ttta<-ttta[!is.na(ttta$Symbol),]
  ttta<-ttta[order(ttta$Symbol),]
  length(ttta$Symbol)
  ttta$y<-rep("1", length(ttta$Symbol)); 
  for (i in 2:length(ttta$Symbol)) { if (ttta$Symbol[i] == ttta$Symbol[i-1]) {ttta$y[i]<-0}}
  ttta[1:3,]; 
  ttta<-ttta[ttta$y==1,]; row.names(ttta)<-ttta$Symbol
  head(ttta)
  all <- rbind(all, ttta)
}
all00 <- all[,-c(dim(all)[2])]
colnames(all00)[23] <- "bnGlobalKDP"

all01 = all00[,c(1,19,22,23,2:18)] 
write.table(all01, file = "DEP_union/Protein_exprCorrected_HvsL_Braak_CDR_CERAD_PLQ_withVotes_uniqueGenSymbol_updated.txt",
             sep="\t", row.names = F, quote = F)

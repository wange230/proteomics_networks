## literature for co-occurrance of a KDP and "Alzheimer's disease"
rm(list=ls())
wd ="C:/Users/wange13/Desktop/Proteomics.manuscript"
setwd(wd)
library(RISmed)
#dir.create("RCG_literature_all")
dir = "RCG_literature_all_KDP"
if (dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}

setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/revision/revision_R1/supplementary/supplementary_tables/supplementary_tables") 

library(readxl)

dataset <- read_excel("Table S4.xlsx", sheet = 9)
dat = data.frame(dataset)
head(dat)
gen = unique(dat$Gene.Symbol)

setwd(wd)
for (gene in gen) {
keywords = paste0(gene, c(" and Alzheimer's disease"))
#keywords = paste0("MSN", c(" and Alzheimer's disease")) # for test
dfff <- NULL
for (kk in keywords) {
res <- EUtilsSummary(kk, type="esearch", db="pubmed", datetype='pdat') 
#,mindate=2000, maxdate=2024, retmax= 200)  # limit of years
#QueryCount(res)
if (QueryCount(res) == 0 ) {next}
tle<-ArticleTitle(EUtilsGet(res))

au <- Author(EUtilsGet(res))
Sys.sleep(10)
abs <- AbstractText(EUtilsGet(res))
Sys.sleep(10)
pmid <- PMID(EUtilsGet(res))
Sys.sleep(10)
y <- YearPubmed(EUtilsGet(res))
Sys.sleep(10)
r <- YearReceived(EUtilsGet(res))

aut <- c();yr  <- c();pmd <- c();Abstract <- list(); Title = list()

for (nn in 1:length(y)) {
aut <- c(aut,au[[nn]]$LastName[1]) 
yr <- c(yr,y[nn])
pmd <- c(pmd,pmid[nn])
Abstract[[nn]] <- abs[nn]
Title[[nn]] = unlist(tle[[nn]])
}

df <- data.frame(Search.words = rep(kk, length(yr)), 
                 Author = aut,
                 Year = yr, 
                 Title = cbind(Title), 
                 Abstract = cbind(Abstract), 
                 PMID = pmd)
head(df)

if(is.null(df)) {next}
if (dim(df)[1] == 1) { df = rbind(df,df)}  
dff <- as.data.frame(apply(df,2,as.character))
dff$Abstract <- gsub('"',"", dff$Abstract)
dff$Abstract <- gsub("", "", gsub("c\\(", "", gsub ("\\)","", dff$Abstract)))
dff$Title <- gsub('"',"", dff$Title)
dff$Title <- gsub("", "", gsub("c\\(", "", gsub ("\\)","", dff$Title)))
dff$Number.literature <- rep(dim(dff)[1],dim(dff)[1])
if(is.null(dff)) {next}
dfff <- rbind(dfff, dff)
rm(aut,yr,pmd,Abstract,Title)
}
# output
#write.table(dfff, file=paste0(paste0("RCG_literature_all_KDP/Resilience_RCG_AD_", gene),"_literature_12_2024.txt"),sep="\t", row.names = FALSE, quote = F)
print(kk)
}
# over
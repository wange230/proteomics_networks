# summary of literature
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/RCG_literature_all_KDP"
setwd(wd)
library(dplyr)
fls = list.files()
head(fls)

fn = fls[grepl(".txt", fls)]
head(fn,20)
prot = gsub("Resilience_RCG_","",gsub("_literature.*","", fn))
head(prot)

df = NULL
nams = c(); num = c()

for(ii in 1:length(fn)){
#  ii = 1
if(file.size(fn[ii]) > 2) {
dat = read.table(fn[ii], sep="\t", stringsAsFactors = F, header = T, quote = '"', fill = T)
head(dat)
if (nrow(dat) == 2 & dat$PMID[1] == dat$PMID[2]){num = c(num, 1)} else {dat = dat %>% distinct(PMID,.keep_all = T); num = c(num, nrow(dat));rm(dat)}
} else{num = c(num, 0)}
nams = c(nams, prot[ii])
print(ii)
}
df = data.frame(KDP= nams, num_co_exist_with_AD = num )
write.table(df, file="co_existance_KDP_AD_in_pubmed.tsv", sep="\t", row.names = F,quote = F )

# combining with KDP rankings
setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/revision/revision_R1/supplementary/supplementary_tables/supplementary_tables") 
library(readxl)
dataset <- read_excel("Table S4.xlsx", sheet = 9)
dat = data.frame(dataset)
head(dat)

datt = dat %>% distinct(Gene.Symbol, .keep_all = T)
head(datt)

df$Gene.Symbol = gsub("AD_","",df$KDP)
length(intersect(datt$Gene.Symbol,df$Gene.Symbol ))
setdiff(datt$Gene.Symbol, df$Gene.Symbol)
# "NA"   "TPM3" "CS"   "HDGF" "SPR"  "TSPO" "RPS3"

head(df)
head(datt)
final = merge(datt, df, by.x = "Gene.Symbol", by.y = "Gene.Symbol", all = F)
setwd(wd)
#write.table(final, file="co_existance_KDP_AD_in_pubmed_withKDPranking.tsv", sep="\t", row.names = F,quote = F )

cor.test(final$KDP.ranking.score, final$num_co_exist_with_AD)
cor.test(final$KDP.ranking.Order, final$num_co_exist_with_AD)
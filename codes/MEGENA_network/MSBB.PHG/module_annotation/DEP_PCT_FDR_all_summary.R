##
rm(list=ls())
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/PHG_traits_cor")

corr <- read.table("proteomics.allAdj_cor_traits_PHG_withID.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
corr.fdr <- corr[corr$p.BH < 0.05, c(1,9)]
head(corr.fdr)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019")
fns <- list.files()[grepl("final.tsv",list.files())]

traits = c("Braak", "CDR","CERAD","PlaqueMean")
dff <- NULL
for (nn in 1:length(traits)) {
fns.braak <- fns[grepl(traits[nn], fns)]
braak <- c()
for (fl in fns.braak) {
  dat <- read.table(fl, sep ="\t", stringsAsFactors = FALSE, header = TRUE)
  dat.sel = dat[abs(dat$logFC) > log2(1.1) & dat$adj.P.Val < 0.05,]
  braak <- unique(c(braak,row.names(dat.sel)))
}
length(braak)

df <- data.frame(protein = braak, module = rep(traits[nn], length(braak)))
dff <- rbind(dff, df)
}
head(dff); tail(dff)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new")
dir = "DEP_PCT_all"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
write.table(corr.fdr, file = "DEP_PCT_all/PCT_FDR_summary.txt", sep="\t", row.names = FALSE, quote = F)
write.table(dff, file = "DEP_PCT_all/DEP_FDR_summary.txt", sep="\t", row.names = FALSE, quote = F)






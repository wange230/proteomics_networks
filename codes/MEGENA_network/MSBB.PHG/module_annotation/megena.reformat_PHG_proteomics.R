##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/MEGENA_new"
setwd(wd)

setwd("C:/Users/wange13/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics")
fdata <- read.table("MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj.tsv",
                    header = TRUE, sep="\t", quote ='"', stringsAsFactors = FALSE)
fdata <- fdata[, 1:2]
head(fdata)

setwd(wd)
megena <- read.table("MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_forMEGENA.tsv--MEGENA2WGCNA.tsv",
                     sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(megena)
final <- merge(fdata, megena, by.x = 2, by.y = 1, all = FALSE)
head(final)
reformat <- final[,c(1:2,8:9)]
head(reformat)
reformat$module <- gsub(".*_","M",reformat$module)
reformat$is.hub[is.na(reformat$is.hub)] <- ""
head(reformat)

setwd(wd)
dir = "MEGENA_reformat_protein"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}
write.table(reformat, file="MEGENA_reformat_protein/MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat.tsv", 
            sep="\t", quote = FALSE, row.names = FALSE)

mod <- c(); parent <- c(); num <- c()
splt <- split(final, final$module)
#splt[[1]][1,5]; splt[[1]][1,9]
for (ii in 1:length(splt)) {
  mod <- c(mod, splt[[ii]][1,9])
  parent <- c(parent, splt[[ii]][1,5])
  num <- c(num, dim(splt[[ii]])[1])
}
df <- data.frame(module = mod, mod.parent = parent, module.size = num)
df$module <- gsub(".*_","M",df$module)
df$mod.parent <- gsub(".*_","M",df$mod.parent)
head(df)

write.table(df, file="MEGENA_reformat_protein/MEGENA_MSBB.PHG.proteomics.PMI_Age_race_sex_batch_adj_MEGENA2WGCNA_reformat_module_summary.tsv", 
            sep="\t", quote = FALSE, row.names = FALSE)


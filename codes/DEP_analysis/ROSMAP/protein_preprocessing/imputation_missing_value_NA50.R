##
rm(list=ls())
library(impute)
setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/output")
tbl <- read.csv("C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400.csv")
dim(is.na(tbl))
row.names(tbl) <-tbl$X
tbl00 <- tbl[,-1]
drop <- rowSums(is.na(tbl00)) >= 0.5*dim(tbl)[2]
drop[1:5]
mat <- as.matrix(tbl00[!drop,])
mat[1:5,1:5]; dim(mat)
exp.25 <- impute.knn(mat, k=50)
df.25 <- data.frame(exp.25$data)

dir = "imputed"
if(dir.exists(dir)){print("This dir has already existed!!!")} else{dir.create(dir)}
write.table(df.25,file="imputed/C2.median_polish_corrected_log2(abundanceRatioCenteredOnMedianOfBatchMediansPerProtein)-8817x400_50na.txt",
            sep="\t", quote = F)

##
rm(list=ls())
setwd("C:/Users/wange13/Desktop/mouse5fad")

mon = c(3,6,12)
dff <- NULL
for (mm in 1:length(mon)) {
deg <- read.table(paste0(paste0("mouse_FAD_vs_WT_DEP/mouse_FAD_vs_WT_month", mon[mm]), "_Sex_adj_DEP.txt"),
                  sep = "\t", stringsAsFactors = FALSE, header = TRUE)
head(deg)

deg.sel.pos <- deg$Protein.Accession[deg$adj.P.Val < 0.05&deg$logFC>0]
head(deg.sel.pos)

deg.sel.neg <- deg$Protein.Accession[deg$adj.P.Val < 0.05&deg$logFC<0]
head(deg.sel.neg)

df <- data.frame(gene.id = c(deg.sel.neg,deg.sel.pos), 
                 module = c(rep(paste0(paste0("month", mon[mm]),".neg"),length(deg.sel.neg)),
                            rep(paste0(paste0("month",mon[mm]),".pos"),length(deg.sel.pos))))
head(df)
dff <- rbind(dff,df)

}
head(dff)
getwd()
write.table(dff, file="summary_DEP_Sex_adj_over_age.txt", row.names = FALSE, sep="\t")


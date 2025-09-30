##
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019/count_DEP_FC11_FDR005_withSymbol/hyperGeTest"
setwd(wd)
library(stringr)
library(ggplot2)
library(reshape2)
dir = "plot"
if(dir.exists(dir)) {print("This dir has already existed!!!")} else {dir.create(dir)}

go <- read.table("EnrichmentTable_PHG.proteomics_cellType500topHuman.txt",
                 sep="\t", stringsAsFactors = F, header = T)
head(go)
go$Set1 <- gsub("PlaqueMean","PLQ",gsub("negative","Dn",gsub("positive","Up",go$Set1)))
go.sel <- go[,c(1,2,9)]
head(go.sel)

df <- cbind(unique(go.sel$Set2))
splt <- split(go.sel, go.sel$Set1)
names(splt)

for (nn in 1:length(splt)) {
dat <- splt[[nn]][,c(2,3)]
colnames(dat)[2] <- names(splt)[nn]
df <- merge(df, dat, by.x= 1, by.y =1, all.x =T)
}

df[is.na(df)] <- 1
df.row.names <- df$V1
row.names(df) <- df$V1
dff <- df[,-1]
mat <- -log10(as.matrix(dff))
dff <- data.frame(mat)
head(dff)
row.names(dff) = str_to_title(row.names(dff))
#head(mtcars)
datt <- melt(dff)
datt$cell.Type <- rep(row.names(dff), dim(go.sel)[1]/6)
head(datt)
datt = datt[!datt$cell.Type%in%"Opc",]
table(datt$variable)
datt$variable = gsub("HvsL_Braak.Dn","High vs. Low.Down",datt$variable)
datt$variable = gsub("HvsL_Braak.Up","High vs. Low.Up",datt$variable)
datt$variable = gsub("HvsL_CDR.Dn","Demented vs. Nondemented.Down",datt$variable)
datt$variable = gsub("HvsL_CDR.Up","Demented vs. Nondemented.Up",datt$variable)
datt$variable = gsub("HvsL_CERAD.Dn","DefiniteAD vs. NL.Down",datt$variable)
datt$variable = gsub("HvsL_CERAD.Up","DefiniteAD vs. NL.Up",datt$variable)
datt$variable = gsub("HvsL_PLQ.Dn","Severe vs. Normal.Down",datt$variable)
datt$variable = gsub("HvsL_PLQ.Up","Severe vs. Normal.Up",datt$variable)
datt$variable = gsub("HvsM_Braak.Dn","High vs. Medium.Down",datt$variable)
datt$variable = gsub("HvsM_Braak.Up","High vs. Medium.Up",datt$variable)
datt$variable = gsub("HvsM_CDR.Dn","Demented vs. Impaired.Down",datt$variable)
datt$variable = gsub("HvsM_CDR.Up","Demented vs. Impaired.Up",datt$variable)
datt$variable = gsub("HvsM_CERAD.Dn","DefiniteAD vs. IndefiniteAD.Down",datt$variable)
datt$variable = gsub("HvsM_CERAD.Up","DefiniteAD vs. IndefiniteAD.Up",datt$variable)
datt$variable = gsub("HvsM_PLQ.Dn","Severe vs. Mild.Down",datt$variable)
datt$variable = gsub("HvsM_PLQ.Up","Severe vs. Mild.Up",datt$variable)
datt$variable = gsub("MvsL_CERAD.Up","IndefintedAD vs. NL.Up",datt$variable)
datt$variable = gsub("MvsL_PLQ.Up","Mild vs. Normal.Up",datt$variable)

datt$variable = factor(datt$variable, 
                       levels = c("High vs. Low.Down",
                                  "Demented vs. Nondemented.Down",
                                  "DefiniteAD vs. NL.Down",
                                  "Severe vs. Normal.Down",
                                  "High vs. Medium.Down",
                                  "Demented vs. Impaired.Down",
                                  "DefiniteAD vs. IndefiniteAD.Down",
                                  "Severe vs. Mild.Down",
                                  "High vs. Low.Up",
                                  "Demented vs. Nondemented.Up",
                                  "DefiniteAD vs. NL.Up",
                                  "Severe vs. Normal.Up",
                                  "High vs. Medium.Up",
                                  "Demented vs. Impaired.Up",
                                  "DefiniteAD vs. IndefiniteAD.Up",
                                  "Severe vs. Mild.Up",
                                  "IndefintedAD vs. NL.Up",
                                  "Mild vs. Normal.Up")
)
colnames(datt)[2] = "-log10(p.adj)"

pdf(file = "plot/PHG.prots.DEP.enrichemnt_over_celltype_06_17_2024.pdf",
      width = 11.0, height = 4.5)
ggplot(datt, aes(variable, cell.Type)) +
  geom_tile(aes(fill = `-log10(p.adj)`), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "", y = "")  +
  theme(axis.text.x = element_text(colour="grey20",size=9.0,
                                   angle= 45,hjust=0.9,vjust=0.95,face="plain"), #original 2
        axis.text.y = element_text(colour="grey20",size=9.0,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=6,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=6,angle=90,hjust=.5,vjust=.5,face="plain"), #original 4
        legend.position="right",
        legend.title = element_text(color = "blue", size = 15),
        legend.text = element_text(color = "red", size = 10),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x =element_line(size=0.5),
        axis.ticks.y =element_line(size=0.5),
        axis.ticks.length=unit(0.01,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black") )
dev.off()
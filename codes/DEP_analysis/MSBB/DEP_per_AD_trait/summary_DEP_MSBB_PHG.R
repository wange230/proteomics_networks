#
rm(list=ls())
wd = "C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019/count_DEP_FC11_FDR005"
setwd(wd)
library(tidyverse)
library(ggplot2)
library(reshape2)

dat = read.table("DEProtein_final_11_12_2019_summary.txt", 
                 sep="\t", stringsAsFactors = F, header = T)
head(dat)
row.names(dat) = dat$comparison
dat = dat[,-ncol(dat)]
head(dat)
datt <- melt(dat)
head(datt)
datt$sign = rep(row.names(dat), 12)
head(datt)
colnames(datt) = c("Comparisons","Numbers","Changes")
df = datt
head(df)
df$Comparisons = gsub("HvsL","ADvsNL",gsub("HvsM","ADvsMCI",gsub("MvsL","MCIvsNL", df$Comparisons)))
df$Comparisons = gsub("_CDR","(CDR)",gsub("_CERAD","CERAD",gsub("_PlaqueMean","(PlaqueMean)",gsub("_Braak","(Braak)",df$Comparisons))))

df$Comparisons = gsub("ADvsMCI\\(Braak\\)","High vs. Medium", df$Comparisons)
df$Comparisons = gsub("ADvsMCI\\(CDR\\)","Demented vs. Impaired", df$Comparisons)
df$Comparisons = gsub("ADvsMCI\\(PlaqueMean\\)","Severe vs. Mild", df$Comparisons)
df$Comparisons = gsub("ADvsMCICERAD","DefiniteAD vs. IndefiniteAD", df$Comparisons)
df$Comparisons = gsub("ADvsNL\\(Braak\\)","High vs. Low", df$Comparisons)
df$Comparisons = gsub("ADvsNL\\(CDR\\)","Demented vs. Nondemented", df$Comparisons)
df$Comparisons = gsub("ADvsNL\\(PlaqueMean\\)","Severe vs. Normal", df$Comparisons)
df$Comparisons = gsub("ADvsNLCERAD","DefiniteAD vs. NL", df$Comparisons)
df$Comparisons = gsub("MCIvsNL\\(Braak\\)","Medium vs. Low", df$Comparisons)
df$Comparisons = gsub("MCIvsNL\\(CDR\\)","Impaired vs. Nondemented", df$Comparisons)
df$Comparisons = gsub("MCIvsNL\\(PlaqueMean\\)","Mild vs. Normal", df$Comparisons)
df$Comparisons = gsub("MCIvsNLCERAD","IndefiteAD vs. NL", df$Comparisons)
table(df$Comparisons)
df$Comparisons = factor(df$Comparisons, levels=c("Medium vs. Low",  
                                                 "High vs. Medium",
                                                 "High vs. Low" , 
                                                 "Impaired vs. Nondemented",
                                                 "Demented vs. Impaired", 
                                                 "Demented vs. Nondemented",
                                                 "IndefiteAD vs. NL",
                                                 "DefiniteAD vs. IndefiniteAD",
                                                 "DefiniteAD vs. NL",
                                                 "Mild vs. Normal",
                                                 "Severe vs. Mild",
                                                 "Severe vs. Normal" 
  
))

df$Changes = gsub("Total.no","Sum",gsub("Positive.no","Up",gsub("Negative.no","Down", df$Changes)))
head(df)
df$Changes = factor(df$Changes, levels=c("Down","Up","Sum"))
g = ggplot(df, aes(Comparisons,Numbers, fill=Changes)) + 
   geom_bar(stat="identity", position=position_dodge())
g1 = g + 
labs(
  x = "",
  y = "Numbers of DEPs",
  ) +
  scale_fill_manual(
    values = c("blue", "red", "turquoise"))
g1
g2 = g1 +
  theme(axis.text.x = element_text(colour="grey20",size=7.5,angle= 30,hjust=0.9,vjust=1.0,face="bold"), #original 2
        axis.text.y = element_text(colour="grey20",size=10,angle=0,hjust=1,vjust=0,face="bold"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="bold"), #orignal 4
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="bold"),
        legend.position = c(0.5, 0.95),
        #legend.justification = c(0, 1),
        #legend.position="bottom",
        legend.direction="horizontal",
        legend.title =	element_text(size=0),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x =element_line(size=0.35),
        axis.ticks.y =element_line(size=0.35),
        axis.ticks.length=unit(0.1,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black") )
g2

setwd(wd)
pdf(file="Summary_DEP_MSBB_PHG_06_05.pdf", height = 4.5, width = 5.0)
print(g2)
dev.off()
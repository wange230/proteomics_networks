##
rm(list=ls())
wd = "C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP"
dir.create("DEP_DEG_diff")
setwd(wd)
setwd("C:/Users/wange13/Desktop/ROSMAP_meth/DEP")
deg <- read.table("ROSMAP_consensus_DEG.txt",
                  sep="\t", stringsAsFactors = F, header = T)
head(deg)
length(deg$Symbol[deg$comp_DEG == "Up"])+length(deg$Symbol[deg$comp_DEG == "Dn"])
setwd("C:/Users/wange13/Desktop/ROSMAP.proteomics/DEP")
dep <- read.table("ROSMAP_consensus_DEP.txt",
                  sep="\t", stringsAsFactors = F, header = T)
head(dep)
table(dep$symbol%in%deg$Symbol)

#superexactTest
background.signature = intersect(deg$Symbol,dep$symbol)
length(background.signature)

dep.sel <- dep[dep$symbol%in%background.signature,]

deg.sel <- deg[deg$Symbol%in%background.signature,]

#superExactTest
a=as.character(dep.sel$symbol[dep.sel$comp_DEP == "Up"]) 
b=as.character(dep.sel$symbol[dep.sel$comp_DEP == "Dn"])
c=as.character(dep.sel$symbol[dep.sel$comp_DEP == "Neutral"])
d=as.character(deg.sel$Symbol[deg.sel$comp_DEG == "Up"]) 
e=as.character(deg.sel$Symbol[deg.sel$comp_DEG == "Dn"])
f=as.character(deg.sel$Symbol[deg.sel$comp_DEG == "Neutral"])

library(ggVennDiagram)
library(ggplot2)
library("ggsci")
x <- list(`DEP(Up)` = a, `DEG(Up)` = d)
#tiff(file="DEP_DEG_diff/ROSMAP_DEP_upvsDEG_up.tiff",
    # width = 1500, height = 1500, res=300)
p= ggVennDiagram(x, label_color = "red",
              set_size = 7.5,
              label_size = 6.5,
              set_color = "blue",
              label_txtWidth = 80,
              edge_size = 2.5)
p + scale_fill_distiller(palette = "Pastel1")
#dev.off()

x <- list(`DEP(Down)` = b, `DEG(Down)` = e)
#tiff(file="DEP_DEG_diff/ROSMAP_DEP_dnvsDEG_dn.tiff",
 #    width = 1500, height = 1500, res=300)
p= ggVennDiagram(x, label_color = "red",
                 set_size = 7.5,
                 label_size = 6.5,
                 set_color = "blue",
                 label_txtWidth = 80,
                 edge_size = 2.5)
p + scale_fill_distiller(palette = "Pastel1")
#dev.off()

# next superexact


table(a%in%c)
set00 <- list("DEP(up)"=a,"DEP(down)"=b,"DEP(neutral)"=c, "DEG(up)"=d, "DEG(down)"=e,"DEG(neutral)"=f)
class(set00); set00[1:2]

getwd()
#Load the SuperExactTest package.
library(SuperExactTest)
str(set00)
#Perform the super exact test.  
Result1 = supertest(set00, n = length(background.signature))
#Visualize the result in a circular layout.
plot(Result1, degree = 2:6, sort.by = 'size',color.scale.pos="topright", 
     legend.col = 1, keep.empty.intersections=FALSE, maxMinusLog10PValue=20)

#tiff("superExactTest_DEP_DEG_PFC_updated.tiff", height = 2000,width = 2000, res=300 )
plot(Result1, degree = 2:6, sort.by = 'size',color.scale.pos="topright", 
     legend.col = 1, keep.empty.intersections=FALSE, maxMinusLog10PValue=100)
#dev.off()
# plot(Result1, degree = 2:7, sort.by = 'size', legend.col = 1, ylim =c(0,4))

#Visualize the result in a matrix layout.
#tiff("DEP_DEG_diff/superExactTest_DEP_DEG_PFC_updated_landscape_2degree.tiff",
 #    height = 2750, width = 3000, res = 150 )
pdf("DEP_DEG_diff/superExactTest_DEP_DEG_PFC_updated_landscape_2degree.pdf",
     height = 7.5, width = 7.5 )
plot(Result1, Layout = 'landscape', 
     degree = 2,
     sort.by = 'size',
     overlap.size.cex = 1.5,
     ylim =c(0,5500),
     gap.between.track = 3.0,
     cex = 1.5,
     cex.lab = 0.5,
     legend.text.cex = 1.5,
     color.scale.cex = 1.5,
     min.intersection.size=1,
     margin=c(2.5,7.5,5.0,3.5),
     ylab = "")
dev.off()

# Tabulate the analysis result into a file:  
#write.csv(summary(Result1)$Table, file = 'DEP_DEG_diff/PFC_DEG_DEP_overlap.csv', row.names = FALSE)

df <- summary(Result1)$Table
df.sel <- df[df$Degree == 2&!df$Observed.Overlap==0,]
dat <- NULL
for (nn in 1:dim(df.sel)[1]) {
strr = (unlist(strsplit(df.sel[nn,7], split=", ")))
dff <- data.frame(GN =strr , module = rep(df.sel[nn,1], length(strr)))
dat <- rbind(dat,dff)
}
head(dat)
table(dat$module)
intersect(as.character(dat$GN[dat$module == "DEP(neutral) & DEG(neutral)"]),as.character(dat$GN[dat$module == "DEP(up) & DEG(up)"]))
dim(dat)
setwd(wd)
#write.table(dat, file = 'DEP_DEG_diff/PFC_DEG_DEP_combination_signatures.tsv', row.names = FALSE, quote = F, sep = "\t")
# 
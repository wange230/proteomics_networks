#
rm(list=ls())
wd = "C:/Users/wange13/Desktop/Proteomics.manuscript/RCG_literature_all_KDP"
setwd(wd)
library(ggplot2)
library(ggpubr)
df = read.table("co_existance_KDP_AD_in_pubmed_withKDPranking.tsv",sep="\t",stringsAsFactors = F, header = T)

cor.test(df$KDP.ranking.score[df$KDP.ranking.score>0], df$num_co_exist_with_AD[df$KDP.ranking.score>0])
cor.test(df$KDP.ranking.Order[df$KDP.ranking.score>0], df$num_co_exist_with_AD[df$KDP.ranking.score>0])

hist(df$num_co_exist_with_AD)

df$num_co_exist_with_AD_categories = ifelse(df$num_co_exist_with_AD == 0, "0", ifelse(df$num_co_exist_with_AD<6, "5", ifelse(df$num_co_exist_with_AD<11, "10", ifelse(df$num_co_exist_with_AD<21, "20", ifelse(df$num_co_exist_with_AD<51,"50", ifelse(df$num_co_exist_with_AD<100, "100", ifelse(df$num_co_exist_with_AD<201, "200",">200")))))))
table(df$num_co_exist_with_AD_categories)
dat = as.data.frame(table(df$num_co_exist_with_AD_categories))
head(dat)
colnames(dat) = c("Categories", "Counts")
dat$Categories = factor(dat$Categories, levels = c(0,5,10,20,50,100,200, ">200"))

p1 = ggplot(dat, aes(x=Categories,y=Counts)) +
     geom_bar(stat = "identity",fill = "#00AFBB") +
     #ggtitle("Number of papers") +
     coord_flip() +
     xlab("Co-occurrance Categories") + ylab("Counts of KDPs")+
     theme_linedraw()
p1
#pdf(file="viz_counts_paper.pdf", width = 4.5, height = 3.5)
#print(p1)
#dev.off()

# not run code below, just for data exploration
df$num_co_exist_with_AD_categories = factor(df$num_co_exist_with_AD_categories, levels = c(0,5,10,20,50,100,200, ">200"))

g = ggplot(df, aes(x = num_co_exist_with_AD)) +
    geom_histogram(color="blue",
                   fill="lightblue",
                   binwidth = 20,
                   #bins = 10,
                   boundary = 0) +
    ggtitle("Frequency histogram")
g


g1 = ggplot(df, aes(x = num_co_exist_with_AD)) +
     geom_histogram(aes(y = after_stat(density)),
                    color = "blue",
                    fill = "lightblue", binwidth = 200) +
  geom_function(fun = dnorm, col = "red", args = list(mean(df$num_co_exist_with_AD), sd(df$num_co_exist_with_AD)), lwd = 1) +
  labs(title = "Normal curve overlaid")
g1

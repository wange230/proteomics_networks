##
rm(list=ls())
setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019/count_DEP_FC11_FDR005_withSymbol/MSigDB")
library(ggplot2)
library(reshape2)
library(stringr)
go <- read.table("DEP_final_02_21_2020_summary_OntologyTop10-sortByModule.xls",
                 sep="\t", stringsAsFactors = F, header = T)
head(go)
go$module <- gsub("PlaqueMean","PLQ",gsub("negative","Dn",gsub("positive","Up",go$module)))

go.sel <- go[!go$module%in%c("MvsL_CERAD.Up","MvsL_PLQ.Up"),c(1,4,10)]
head(go.sel)

# only use HvsL
go.sel = go.sel[!grepl("HvsM_",go.sel$module),]


splt00 <- split(go.sel, go.sel$module)
tbl <- NULL
for (ww in 1:length(splt00)) {
#splt00[[1]]
tbl <- rbind(tbl,splt00[[ww]][1:10,])
}

go.sel <- tbl
go.sig <- unique(go.sel$Gene.Category)
#go.sig[52] <- "RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS"
df <- cbind(go.sig)

splt <- split(go.sel, go.sel$module)
names(splt)

for (nn in 1:length(splt)) {
dat <- splt[[nn]][,c(2,3)]
colnames(dat)[2] <- names(splt)[nn]
df <- merge(df, dat, by.x= 1, by.y =1, all.x =T)
}

df[is.na(df)] <- 1
df.row.names <- df$go.sig
row.names(df) <- df$go.sig
dff <- df[,-1]
mat <- -log10(as.matrix(dff))
dff <- data.frame(mat)
dff <- dff[order(row.names(dff)),]
row.names(dff) <- gsub("GO_","",row.names(dff))
head(dff)
row.names(dff) = str_to_title(row.names(dff))
row.names(dff) = gsub("_"," ",row.names(dff))
head(dff)
row.names(dff)[20]<- "Mitochondrial respiratory chain complex I"
row.names(dff)[21]<- "NADH dehydrogenase complex"
#row.names(dff)[54]<- "KEGG ECM receptor interaction"
#row.names(dff)[55]<- "KEGG parkinsons disease"
row.names(dff)[37]<- "KEGG ribosome"
row.names(dff)[39]<- "TCA cycle & respiratory electron transport"
#row.names(dff)[59]<- "Respiratory electron transport ATP synthesis"
#row.names(dff)[60]<- "YGTCCTTGR unknown"

datt <- melt(dff)
datt$go <- rep(row.names(dff), dim(go.sel)[1]/10)
head(datt)
datt$value[datt$value > 30] = 30
colnames(datt)[2] = "-log10(p.adj)"
datt = datt[!datt$go%in%"Respiratory electron transport ATP synthesis",]

datt$variable = gsub("HvsL_Braak.Dn","High vs. Low.Down",datt$variable)
datt$variable = gsub("HvsL_Braak.Up","High vs. Low.Up",datt$variable)
datt$variable = gsub("HvsL_CDR.Dn","Demented vs. Nondemented.Down",datt$variable)
datt$variable = gsub("HvsL_CDR.Up","Demented vs. Nondemented.Up",datt$variable)
datt$variable = gsub("HvsL_CERAD.Dn","DefiniteAD vs. NL.Down",datt$variable)
datt$variable = gsub("HvsL_CERAD.Up","DefiniteAD vs. NL.Up",datt$variable)
datt$variable = gsub("HvsL_PLQ.Dn","Severe vs. Normal.Down",datt$variable)
datt$variable = gsub("HvsL_PLQ.Up","Severe vs. Normal.Up",datt$variable)
#datt$variable = gsub("HvsM_Braak.Dn","High vs. Medium.Down",datt$variable)
#datt$variable = gsub("HvsM_Braak.Up","High vs. Medium.Up",datt$variable)
#datt$variable = gsub("HvsM_CDR.Dn","Demented vs. Impaired.Down",datt$variable)
#datt$variable = gsub("HvsM_CDR.Up","Demented vs. Impaired.Up",datt$variable)
#datt$variable = gsub("HvsM_CERAD.Dn","DefiniteAD vs. IndefiniteAD.Down",datt$variable)
#datt$variable = gsub("HvsM_CERAD.Up","DefiniteAD vs. IndefiniteAD.Up",datt$variable)
#datt$variable = gsub("HvsM_PLQ.Dn","Severe vs. Mild.Down",datt$variable)
#datt$variable = gsub("HvsM_PLQ.Up","Severe vs. Mild.Up",datt$variable)

datt$variable = factor(datt$variable, 
                       levels = c("High vs. Low.Down",
                                  "Demented vs. Nondemented.Down",
                                  "DefiniteAD vs. NL.Down",
                                  "Severe vs. Normal.Down",
                                  #"High vs. Medium.Down",
                                  #"Demented vs. Impaired.Down",
                                  #"DefiniteAD vs. IndefiniteAD.Down",
                                  #"Severe vs. Mild.Down",
                                  "High vs. Low.Up",
                                  "Demented vs. Nondemented.Up",
                                  "DefiniteAD vs. NL.Up",
                                  "Severe vs. Normal.Up"#,
                                  #"High vs. Medium.Up",
                                  #"Demented vs. Impaired.Up",
                                  #"DefiniteAD vs. IndefiniteAD.Up",
                                  #"Severe vs. Mild.Up"
                                  )
                       )
old.lels = c("Cell projection part" 
             , "Cell projection"
             , "Cell junction"
             , "Neuron part"                                 
             , "Neuron projection"
             , "Axon"
             , "Synapse"                                     
             , "Synapse part"                                
             , "Synaptic signaling"
             , "Postsynapse"                                 
             , "Presynapse"
             , "Somatodendritic compartment"  
             , "Ribosomal subunit"
             , "Cytoplasmic vesicle part"
             , "Exocytic vesicle membrane"
             , "Respiratory chain"
             , "Cytosolic large ribosomal subunit"
             , "NADH dehydrogenase complex" 
             , "Mitochondrial respiratory chain complex I"  
             , "Oxidative phosphorylation"                   
             , "Oxidoreductase complex" 
             , "Membrane protein complex"
             , "Reactome respiratory electron transport" 
             , "Electron transport chain"
             , "Cellular respiration" 
             , "KEGG parkinsons disease"                     
             , "KEGG ribosome"
             , "Regulation of exocytosis"   
             , "YGTCCTTGR unknown"
             , "TCA cycle & respiratory electron transport" 
             
             
             , "Actin binding"                               
             , "Actin cytoskeleton"                          
             , "Actin filament based process"                
             , "Anchoring junction"                                                                     
             , "Basement membrane"                           
             , "Biological adhesion"       
             , "Cell substrate junction"                     
             , "Cellular response to organic substance"      
             , "Cytoskeletal protein binding"                
             , "Defense response"                            
             , "Regulation of cell adhesion"                   
             , "Extracellular matrix"                        
             , "Extracellular matrix component"              
             , "Extracellular structure organization"        
             , "Identical protein binding"                   
             , "Immune effector process"                     
             , "Immune response"                             
             , "Immune system process"                       
             , "Innate immune response" 
             , "Positive regulation of immune system process"
             , "KEGG ECM receptor interaction"    
             , "Oxidoreductase activity"
             , "Protein complex binding"                     
             , "Proteinaceous extracellular matrix"          
             , "S100 protein binding"                        
             , "Sarcolemma"                                  
             , "Single organism catabolic process"           
             , "Single organism cell adhesion"               
             , "Small molecule metabolic process")
setdiff(old.lels, datt$go)
datt$go = factor(datt$go, 
                 levels = c(#"Cell projection part", 
                             "Cell projection"
                            , "Cell junction"
                            , "Neuron part"                                 
                            , "Neuron projection"
                            #, "Axon"
                            , "Synapse"                                     
                            , "Synapse part"                                
                            , "Synaptic signaling"
                            #, "Postsynapse"                                 
                            #, "Presynapse"
                            #, "Somatodendritic compartment"  
                            , "Ribosomal subunit"
                            #, "Cytoplasmic vesicle part"
                            #, "Exocytic vesicle membrane"
                            , "Respiratory chain"
                            , "Cytosolic large ribosomal subunit"
                            , "NADH dehydrogenase complex" 
                            , "Mitochondrial respiratory chain complex I"  
                            , "Oxidative phosphorylation"                   
                            , "Oxidoreductase complex" 
                            , "Membrane protein complex"
                            , "Reactome respiratory electron transport" 
                            , "Electron transport chain"
                            , "Cellular respiration" 
                            #, "KEGG parkinsons disease"                     
                            , "KEGG ribosome"
                            #, "Regulation of exocytosis"   
                            #, "YGTCCTTGR unknown"
                            , "TCA cycle & respiratory electron transport" 
                            
                              
                            , "Actin binding"                               
                            #, "Actin cytoskeleton"                          
                            , "Actin filament based process"                
                            , "Anchoring junction"                                                                     
                            , "Basement membrane"                           
                            , "Biological adhesion"       
                            , "Cell substrate junction"                     
                            , "Cellular response to organic substance"      
                            , "Cytoskeletal protein binding"                
                            #, "Defense response"                            
                            #, "Regulation of cell adhesion"                   
                            , "Extracellular matrix"                        
                            , "Extracellular matrix component"              
                            , "Extracellular structure organization"        
                            , "Identical protein binding"                   
                            #, "Immune effector process"                     
                            #, "Immune response"                             
                            , "Immune system process"                       
                            #, "Innate immune response" 
                            #, "Positive regulation of immune system process"
                            #, "KEGG ECM receptor interaction"    
                            , "Oxidoreductase activity"
                            , "Protein complex binding"                     
                            , "Proteinaceous extracellular matrix"          
                            #, "S100 protein binding"                        
                            #, "Sarcolemma"                                  
                            , "Single organism catabolic process"           
                            , "Single organism cell adhesion"               
                            , "Small molecule metabolic process"))



# output
pdf(file = "PHG.proteomics.DEP.GO.top10_06_11_2025.pdf",width = 2.0, height = 4.75)
ggplot(datt, aes(variable, go)) +
  geom_tile(aes(fill = `-log10(p.adj)`), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "", y = "")  + 
  theme(#axis.text.x = element_text(colour="grey20",size=11,
                                  # angle= 30,hjust=1.05,vjust=1.05,face="plain"), #original 2
       # axis.text.y = element_text(colour="grey20",size=11,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=0,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="grey20",size=0,angle=90,hjust=.5,vjust=.5,face="plain"), #original 4
        legend.position="none",
        #legend.title = element_text(color = "blue", size = 10),
        #legend.text = element_text(color = "red", size = 5),
        axis.line.x = element_line(size=0.5),
        axis.line.y = element_line(size=0.5),
        axis.ticks.x =element_line(size=0.5),
        axis.ticks.y =element_line(size=0.5),
        axis.ticks.length=unit(0.0,"cm"),plot.margin = margin(0,0,0,0,"pt"),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black") )
dev.off()
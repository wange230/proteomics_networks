##
rm(list=ls())
setwd("C:/Users/wange13/Desktop/mouse5fad")
# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  #humanx <- unique(genesV2[, 2])
  humanx <- genesV2[, ]
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  #humanx <- unique(genesV2[, ])
  humanx <- unique(genesV2[, ])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

#genes <- convertMouseGeneList(humGenes)
#musGenes <- c("Hmmr", "Tlx3", "Cpeb4")
#mo2hu.homologue <- convertMouseGeneList(musGenes)
#hu2mo.homologue <- convertHumanGeneList(c("TLX3" , "HMMR" , "CPEB4"))

setwd("C:/Users/wange13/Desktop/mouse5fad")
dat <- read.csv("5xFAD_mouse_transcriptome_PengLab_at_StJude_v1.1.0.csv")
musGenes=as.character(dat$GN)
mo2hu.homologue <- convertMouseGeneList(musGenes)
dim(mo2hu.homologue); head(mo2hu.homologue)

write.table(mo2hu.homologue, file= "Mouse2Human_gene_symbal_conversion_00.txt",
            sep="\t", row.names = FALSE)


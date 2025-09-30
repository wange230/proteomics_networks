#
rm(list=ls())
wd = "/sc/arion/projects/zhangb03a/shared/DET_WINA_MEGENA/MSBB_protein_to_ROSMAP_mRNA_preservation/"
setwd(wd)

setwd("/sc/arion/projects/zhangb03a/shared/DET_WINA_MEGENA/PHG_ROSMAP_proteinMEGENA_module_conservation/MSBB_ROSMAP_conservation/")
load("RSOMAP_MSBB_proteomics_conservation.RData")
rosmap = ros01
dim(rosmap)
rosmap[1:5,1:5]
msbb = msbb01
dim(msbb)
msbb[1:5,1:5]

load("ROSMAP_mRNA_proteomics/ROSMAMP_proteomics_mRNA_data4MEGENA_preservation_analysis.RData")

rosmap.rna = mRNA.megena
dim(rosmap.rna)
rosmap.rna[1:5,1:5]

setwd("/sc/arion/projects/zhangb03a/shared/DET_WINA_MEGENA/PHG_MEGENA_module_conservation/for_preservation/") 
msbb.prot = read.table("PHG_proteomics.allAdj_9209Value_with_uniqueID.txt",
                       sep="\t",stringsAsFactors = F,header = T)
msbb.prot[1:5,1:5]

common.gene = intersect(row.names(msbb.prot), row.names(rosmap.rna))
length(common.gene)
length(unique(common.gene))

msbb.exp.prot = msbb.prot[row.names(msbb.prot)%in%common.gene,]
length(row.names(msbb.exp.prot))
length(unique(row.names(msbb.exp.prot)))
rosmap.exp.rna = rosmap.rna[row.names(rosmap.rna)%in%common.gene,]
length(unique(row.names(rosmap.exp.rna)))
length(row.names(rosmap.exp.rna))

setwd(wd)
save(msbb.exp.prot, rosmap.exp.rna, file="data_4_msbb.exp.prot9209_and_rosmap.exp.rna.RData")
# over
# 
rm(list=ls())
setwd("C:/Users/wange13/Desktop/Proteomics.manuscript/proteomics_manuscript/proteomics_manu_v5/adgwas/hyperGeTest")

mm = read.table("EnrichmentTable_PHG_prot_MEGENA_mod50_for_MIT_TWAS.txt",
                  sep="\t",stringsAsFactors = F, header = T)
gen = mm$Overlap.genes[mm$Set1 == "M3"]
gen = unlist(strsplit(gen,";"))

setwd("C:/pc.backup/common_tools_datasets/IGAP_signif_p005")
adgwas = read.table("signature_abeta_adgwas.tsv",
                     sep="\t",stringsAsFactors = F,header = T)
adgwas35 = adgwas$Signature[adgwas$Module == "ADGWAS"] # abeta
abeta = adgwas$Signature[adgwas$Module == "abeta"] # abeta
intersect(gen, adgwas35)
#"APOE"   "PICALM"
intersect(gen, abeta)
# "CAPN2"  "RAB11B" "SNX8"   "ADAM17" "CTNND2" "PTPRC"  "APOE"   "DNAJB6"
#PICALM" "UBQLN1"

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/processed/DEProtein_final_11_12_2019/DEP_union")
dep = read.table("Protein_exprCorrected_HvsL_Braak_CDR_CERAD_PLQ_withVotes_uniqueGenSymbol_updated.txt",
                 sep="\t", stringsAsFactors = F,header = T)
dep00 = unique(dep$Symbol[dep$comp.ADvsNL=="Up"])
dep01 = unique(dep$Symbol[dep$comp.ADvsNL=="Down"])

length(intersect(gen, dep00))
length(intersect(gen, dep01))



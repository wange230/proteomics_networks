rm(list = ls())

make.bsubfile <- function(n.core,job.name,cue.name,project.name,n.hrs,ptile,mem,job.script,mach = "manda")
{
 if (!is.null(mach))
 {
  str.vec <- c("#BSUB -L /bin/bash",paste("#BSUB -n",n.core),paste("#BSUB -R span[ptile=",ptile,"]",sep = ""),paste("#BSUB -R rusage[mem=",mem,"]",sep = ""),paste("#BSUB -J",job.name),paste("#BSUB -o ",job.name,".out",sep = ""),
  paste("#BSUB -e ",job.name,".err",sep = ""),paste("#BSUB -m",mach),paste("#BSUB -q",cue.name),paste("#BSUB -P",project.name),paste("#BSUB -W ",n.hrs,":00",sep = ""),job.script, "\n")
 }else{
  str.vec <- c("#BSUB -L /bin/bash",paste("#BSUB -n",n.core),paste("#BSUB -R span[ptile=",ptile,"]",sep = ""),paste("#BSUB -R rusage[mem=",mem,"]",sep = ""),paste("#BSUB -J",job.name),paste("#BSUB -o ",job.name,".out",sep = ""),
  paste("#BSUB -e ",job.name,".err",sep = ""),paste("#BSUB -q",cue.name),paste("#BSUB -P",project.name),paste("#BSUB -W ",n.hrs,":00",sep = ""),job.script, "\n")
 }
  
 output <- paste(str.vec,collapse = "\n")
 return(output)
}

###############################

root.dir <- "/sc/hydra/projects/zhangb03a/shared/DET_WINA_MEGENA";setwd(root.dir)
data.dir <- "raw_data_PHG_protein"
template.file <- "/sc/hydra/projects/zhangb03a/shared/DET_WINA_MEGENA/tools/R_template_0.txt"
script.dir <- "Scripts_06/"

flag <- c("\\[rootname\\]","\\[datafname\\]","\\[ncore\\]")

MEGENA.ncore = 12
MEGENA.mem = 8000
cue.name <- "premium"
pj.name <- "acc_adineto"
n.hrs <- "138"
####### generate scripts per data file

job.script <- readLines(template.file)
data.files <- list.files(path = data.dir,full.names = T);
# data.files <- data.files[-grep("data\\.annot\\.txt",data.files)]

rfnames <- paste(script.dir,gsub("\\.tsv","\\.R",gsub("^(.*)/","",data.files)),sep = "/")
lsf.files <- paste(script.dir,gsub("\\.tsv","\\.lsf",gsub("^(.*)/","",data.files)),sep = "/")

for (i in 1:length(data.files))
{

 # generate R scripts according to each input file
 new.val <- c(root.dir,data.files[i],MEGENA.ncore)
 file.job <- job.script
 for (j in 1:length(flag))
 {
  file.job <- gsub(flag[j],new.val[j],file.job)
 }

 sink(rfnames[i])
 cat(paste(file.job,collapse = "\n"))
 sink()
 
 R.job <- c(paste("cd ",root.dir,"/",gsub("\\./","",script.dir),sep = ""),"module load R",paste("R CMD BATCH ",gsub("^(.*)/","",rfnames[i]),sep = ""))
 lsf.header <- make.bsubfile(n.core = MEGENA.ncore,job.name = gsub("^(.*)/","",data.files[i]),cue.name = cue.name,project.name = pj.name,n.hrs = n.hrs,ptile = MEGENA.ncore,mem = MEGENA.mem,job.script = R.job,mach = "manda")
 sink(lsf.files[i])
 cat(paste(lsf.header,collapse = "\n"))
 sink()
 
}

bsub.file <- paste(script.dir,"bsub_cmds.txt",sep = "/");
sink(bsub.file)
cat(paste(paste("bsub < ",gsub("^(.*)/","",lsf.files),sep = ""),collapse = "\n"))
sink()

quit(save = "no",status = 0)

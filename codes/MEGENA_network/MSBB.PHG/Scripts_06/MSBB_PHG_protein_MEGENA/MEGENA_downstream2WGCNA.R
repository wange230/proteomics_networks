# convert MEGENA results into WGCNA format
rm(list = ls())
wd <-"/sc/hydra/projects/zhangb03a/shared/DET_WINA_MEGENA/PHG_proteomics_MEGENA"  
setwd(wd)
list.files()
library(MEGENA)

mod.pvalue = hub.pvalue = 0.05
# mod.pvalue = hub.pvalue = 0.01
fns <- list.files()

for (file in fns) {

g <- graph.data.frame(read.delim(file = paste0(file,"/MEGENA_Network.txt"),sep = "\t",
                                 header = TRUE,stringsAsFactors = FALSE),directed = FALSE)
load(paste0(file, "/MEGENA_output.RData"))


module_convert_to_wgcn <- function(MEGENA.output,mod.pval,hub.pval,min.size,max.size)
{     
  summary.output <- MEGENA.ModuleSummary(MEGENA.output,mod.pvalue = mod.pval,hub.pvalue = hub.pval,
                                         min.size = min.size,max.size = max.size,annot.table = NULL,symbol.col = NULL,id.col = NULL,
                                         output.sig = TRUE)
  
  modules <- summary.output$modules
  #####
  is.annted <- any(sapply(modules,function(x) length(grep("\\|",x)) > 0))
  df <- mapply(FUN = function(x,y) data.frame(id = x,module = rep(y,length(x))),x = modules,y = as.list(names(modules)),SIMPLIFY = FALSE)
  df <- do.call('rbind.data.frame',df)
  if (is.annted)
  {
    df <- data.frame(id = df[[1]],gene.symbol = gsub("\\|(.*)","",as.character(df[[1]])),gene.symbol2 = gsub("^(.*)\\|","",as.character(df[[1]])),
                     module = df[[2]])
  }
  
  # update module parent relationship
  if (!is.null(summary.output$module.table))
  {
    modtbl <- summary.output$module.table
    df <- cbind.data.frame(df[,-ncol(df)],module.parent = modtbl$module.parent[match(df$module,modtbl$module.id)],module = df[[ncol(df)]])
    rm(modtbl)
    colnames(df)[1] <- "id"
  }
  
  # update node statistics
  if (!is.null(MEGENA.output$node.summary))
  {
    colnames(df)[1] <- "id"
    # get hub summary
    hubs <- lapply(MEGENA.output$hub.output$module.degreeStat,function(x,hp) as.character(x[[1]])[which(x$pvalue < hp)],hp = hub.pval)
    hvec <- rep(NA,nrow(df))
    for (j in 1:length(hubs)) hvec[which(df[[1]] %in% hubs[[j]] & df$module == names(hubs)[j])] <- "hub"
    df <- cbind.data.frame(df[,-ncol(df)],node.degree = MEGENA.output$node.summary$node.degree[match(df$id,MEGENA.output$node.summary$id)],
                           node.strength = MEGENA.output$node.summary$node.strength[match(df$id,MEGENA.output$node.summary$id)],is.hub = hvec,
                           module = df[[ncol(df)]]);
    rm(MEGENA.output,summary.output)
  }
  rownames(df) <- NULL
  return(df)
}

wgcna <- module_convert_to_wgcn (MEGENA.output=MEGENA.output,mod.pval=mod.pvalue, hub.pval=hub.pvalue, min.size = 10, max.size = vcount(g)/2)

dim(wgcna)

wgcna[1:10,]

wgca <- split(wgcna, wgcna$module); length(wgca); class(wgca); wgca[[18]]

head(wgca)

dir = "MEGENA2WGCNA_10genes"
if(dir.exists("MEGENA2WGCNA_10genes")){print("This dir has already existed!!!")}else {dir.create(dir)}

write.table(wgcna, file = paste0(paste0("MEGENA2WGCNA_10genes/", file), "--MEGENA2WGCNA.tsv"), sep ="\t", row.names = FALSE, quote = FALSE)
}



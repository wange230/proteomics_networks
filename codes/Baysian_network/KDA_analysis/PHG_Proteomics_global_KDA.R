rm(list=ls())
wd <- "C:/Users/wange13/Desktop/MSBB_new_rel"
setwd(wd)

setwd("C:/pc.backup/Desktop/COHORT_data_06_05_2017/MSBB/PHG_Proteomics/bn_network")
ntwk <- read.table("global.BN_cys.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)
ntwk <- ntwk[,1:2]
head(ntwk); dim(ntwk)
fcausalnet <- ntwk

bn <- read.table("global.BN_cys.tsv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
nodes <- union(bn$Source, bn$Target)
finput <- data.frame(id = nodes, module = rep("bn", length(nodes)))
head(finput)
finputlist <- finput

# example 1: key drivers for breast cancer gene modules
setwd(wd)
library(KDA)
#data( breastcausalnet )
#data( breastinputlist )
directed <- TRUE
#directed <- FALSE
layer <- 0
minDsCut <- -1
fgeneinfo <- NULL

# 2. specify the directory for holding analysis results
#
if ( directed )
{
	outputDir <- "KeyDrivers/"
} else
{
	outputDir <- "KeyDrivers-undirected_Proteomics_MEGENA/"
}

dir.create( outputDir )

#
# -----------------------------End of Parameters to be changed --------------------------------------

library( class )
library( cluster )
library( rpart )
#library( sma ) # this is needed for plot.mat below
library( lattice ) # require is design for use inside functions 

memory.size( TRUE )   # check the maximum memory that can be allocated
memory.limit( size = 3800 )   # increase the available memory

################################################################################################
#    1. read in network

cnet <- fcausalnet
cnet <- as.matrix( cnet )
dim( cnet )

totalnodes <- union( cnet[,1] , cnet[,2] )

fname <- "AD"
fname <- paste( fname , "_L" , layer , sep = "" )

################################################################################################
# 2. read in gene lists

listMatrix <- finputlist
dim( listMatrix )
listMatrix <- as.matrix( listMatrix )
listMatrix[1:2,]
ncols <- dim( listMatrix )[2]

modules <- names( table( listMatrix[,ncols] ) )

xkdFall <- paste( outputDir , fname , "_KDx_combined.xls" , sep = "" )
xkdFpara <- paste( outputDir , fname , "_KDx_parameters.xls" , sep = "" )
xkdrMatrix <- NULL
paraMatrix <- NULL

################################################################################################
# 3. process each gene list
#

for ( em in modules )
{
	
#em="green"
	
	print( paste( "*****************" , em , "********************" ) )
	
	esel <- listMatrix[,ncols] == em
	
# remove abnormal gene names
#
	genes <- union( listMatrix[esel,1] , NULL )
	genes <- genes[genes != ""]
	genes <- genes[!is.na( genes )]
	no.genes <- length( genes )
	
	em2 <- replaceString( em , ":" , "" )
	em2 <- replaceString( em2 , " " , "-" )
	
	key2 <- paste( fname , "_KD_" , em2 , sep = "" )
	onetFname <- paste( outputDir , key2 , ".pair" , sep = "" )
	snpFname <- paste( outputDir , key2 , ".snp" , sep = "" )
	kdFname <- paste( outputDir , key2 , "_keydriver.xls" , sep = "" )
	
	if(layer >=1 )
	{
		# expand network by K-hop nearest neighbors layers
		expandNet <- findNLayerNeighborsLinkPairs( linkpairs = cnet , subnetNodes = genes ,
				nlayers = layer , directed = FALSE )
	} 
	else
	{
		# no expansion
		expandNet <- getSubnetworkLinkPairs( linkpairs = cnet , subnetNodes = genes )
	}
	dim( expandNet )
	
	allnodes <- union( expandNet[,1] , expandNet[,2] )
write.table(expandNet, onetFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)
	
	################################################################################################
# 4. keydriver for a given network
#
	
  if (directed)
  {
    ret <- keydriverInSubnetwork( linkpairs = expandNet , signature = genes, background=NULL, directed = directed ,
			         nlayers = 6 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
			         boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }
  else
  {
    ret <- keydriverInSubnetwork( linkpairs = expandNet , signature = genes , directed = directed ,
			         nlayers = 2 , enrichedNodes_percent_cut=-1, FET_pvalue_cut=0.05,
				boost_hubs=T, dynamic_search=T, bonferroni_correction=T, expanded_network_as_signature =F)
  }
	
	if ( is.null( ret ) )
	{
		next
	}
	
	fkd <- ret[[1]]
	parameters <- ret[[2]]
	
	fkd2 <- cbind( rep( em , dim( fkd )[1] ) , fkd )
	xkdrMatrix <- rbind( xkdrMatrix , fkd2 )
	
	nodeDegree = ret[[3]] ##
	nodeDegree[,1] = allnodes[as.integer(nodeDegree[,1])] ###
	
	paraMatrix <- rbind( paraMatrix , c( key2 , parameters ) )
	
	 write.table( fkd , kdFname , sep = "\t" , quote = FALSE , col.names = TRUE , row.names = FALSE )
	
	################################################################################################
# 4. output networks & key drivers for visualization
#
#     Cytoscape output: 1) network file - *_cys.txt 2) node property file: *_cys-nodes.txt
#
	
      nodeprop = configureNodeVisualization(allnodes=allnodes, signature=genes, kdaMatrix=fkd)

      hnList     = nodeprop[[1]] # node subcategpries
      listprop   = nodeprop[[2]] # visual properties for each subcategory
      legend     = nodeprop[[3]] # legend table for visual propertie

      resf = makeSNP(netpairsWtype   = expandNet, 
               edgecolorlevels = c("grey"),
               highlightNodes  = hnList,
               normColor="grey",   highColor=listprop[,1],
               normShape="circle", highShape=listprop[,2],
               normNodeSize ="40",  highNodeSize =listprop[,3],
               normFontSize ="12",  highFontSize =listprop[,4],
               legendtable=legend, snafile=snpFname, kdaMatrix= nodeDegree)
	
} #for (em in modules) {


# save all key drivers
colnames( xkdrMatrix ) <- c( "module" , colnames( fkd ) )
write.table( xkdrMatrix , xkdFall , sep = "\t" , quote = FALSE , col.names = TRUE , 
		row.names = FALSE )     

# save parameters used
#
colnames( paraMatrix ) <- c( "subnet" , colnames( parameters ) )
write.table( paraMatrix , xkdFpara , sep = "\t" , quote = FALSE , col.names = TRUE ,row.names = FALSE )

if( !is.null( fgeneinfo ) )
{
	infoMatrix <- read.delim( fgeneinfo , sep = "\t" , header = TRUE )
	dim( infoMatrix )
	xkdrMatrix2 <- cbind( xkdrMatrix , c( 1:( dim( xkdrMatrix )[1] ) ) )
	
	ic <- dim( infoMatrix )[2] + 1
	merged <- merge( infoMatrix , xkdrMatrix2 , by.x = 1 , by.y = 2 , all.y = TRUE )
	merged <- as.matrix( merged )
	ic2 <- dim( merged )[2]
	xf <- merged[,c( ic , setdiff( 1:ic2 , ic ) )]
	ic3 <- dim( merged )[2]
	mo <- order( as.integer( merged[,ic3] ) )
	 write.table( xf[mo,-ic3] , xkdFall , sep = "\t" , quote = FALSE , col.names = TRUE ,row.names = FALSE )
}

## ------------------------------------- END ------------------------------------------------

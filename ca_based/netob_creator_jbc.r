args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Please supply follwing arguements 1dummydcd 2pdb 3cijresmat 4output_append", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  dcdsmall=args[1]
  pdbf=args[2]
  cijresmat= args[3]
  cutt=0.6
  outapp=args[4]
}
library(bio3d)
library(igraph)

dummydcd<-read.dcd(dcdsmall)
pdb <- read.pdb(pdbf)
ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dummydcd,fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
print("xyzdone")
cij<-dccm(xyz[,ca.inds$xyz])
#from dummy dcd
cmapread=read.table(file=cijresmat,sep=",")
cmapread=as.matrix(sapply(cmapread, as.numeric))
rownames(cmapread)<-cij[,0]
colnames(cmapread)<-cij[0,]

network.amendment <- function(x, membership, collapsem, minus.log=TRUE){

  ## Check for presence of igraph package
  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
     stop("igraph package missing: Please install, see: ?install.packages")
  }

  if(!inherits(x, "cna")){
    stop("Input x must be a 'cna' class object as obtained from cna()")
  }

  if(!is.numeric(membership)){
    stop("Input membership must be a numeric vector")
  }

  if(length(membership) != length(x$communities$membership)){
    stop("Input membership and x$community$membership must be of the same length")
  }
     
  contract.matrix <- function(cij.network, membership,## membership=comms$membership,
                              collapse.method=collapsem, minus.log=TRUE){
    
    ## Function to collapse a NxN matrix to an mxm matrix
    ##  where m is the communities of N. The collapse method
    ##  can be one of the 'collapse.options' below

    ## convert to the original cij values if "-log" was used
    
    if(minus.log){
      cij.network[cij.network>0] <- exp(-cij.network[cij.network>0])
    }
    
    collapse.options=c("max", "median", "mean", "trimmed")
    collapse.method <- match.arg(tolower(collapse.method), collapse.options)

    ## Fill a 'collapse.cij' nxn community by community matrix
    node.num <- max(x$communities$membership)
    if(node.num > 1){
      collapse.cij <- matrix(0, nrow=node.num, ncol=node.num)
      inds <- pairwise(node.num)

      for(i in 1:nrow(inds)) {
        comms.1.inds <- which(membership==inds[i,1])
        comms.2.inds <- which(membership==inds[i,2])
        submatrix <- cij.network[comms.1.inds, comms.2.inds]

        ## Use specified "collapse.method" to define community couplings
        collapse.cij[ inds[i,1], inds[i,2] ] = switch(collapse.method,
                    max = max(submatrix),
                    median = median(submatrix),
                    mean = mean(submatrix),
                    trimmed = mean(submatrix, trim = 0.1))
      }
      
      if(minus.log){
        collapse.cij[collapse.cij>0] <- -log(collapse.cij[collapse.cij>0])
      }
      
      ## Copy values to lower triangle of matrix and set colnames
      collapse.cij[ inds[,c(2,1)] ] = collapse.cij[ inds ]
      colnames(collapse.cij) <- 1:ncol(collapse.cij)
    }
    else{
      warning("There is only one community in the $communities object.
               $community.cij object will be set to 0 in the 
               contract.matrix() function.")

      collapse.cij <- 0
    }

    class(collapse.cij) <- c("dccm", "matrix")
    return(collapse.cij)
  }


  x$communities$membership <- membership

  x$community.cij <- contract.matrix(x$cij, membership, minus.log=minus.log)
  
  cols=vmd_colors()
  
#  if(sum(x$community.cij)>0){
    x$community.network <-  igraph::graph.adjacency(x$community.cij,
                                          mode="undirected",
                                          weighted=TRUE,
                                          diag=FALSE)
        
    ##-- Annotate the two networks with community information
    ## Check for duplicated colors
    if(max(x$communities$membership) > length(unique(cols)) ) {
      warning("The number of communities is larger than the number of unique 
              'colors' provided as input. Colors will be recycled")
    }
  
    ## Set node colors
    igraph::V(x$network)$color <- cols[x$communities$membership]
    igraph::V(x$community.network)$color <- cols[ 1:max(x$communities$membership)]
  
    ## Set node sizes
    igraph::V(x$network)$size <- 1
    igraph::V(x$community.network)$size <- table(x$communities$membership)

#  }

  return(x)
}

mod.select <- function(x,collapsem, thres=0.1) {
   remodel <- community.tree(x, rescale = TRUE)
   n.max = length(unique(x$communities$membership))
   ind.max = which(remodel$num.of.comms == n.max)
   v = remodel$modularity[length(remodel$modularity):ind.max]
   v = rev(diff(v))
   fa = which(v>=thres)[1] - 1
   ncomm = ifelse(is.na(fa), min(remodel$num.of.comms), n.max - fa)
   print (ncomm)
   ind <- which(remodel$num.of.comms == ncomm)
   network.amendment(x, remodel$tree[ind, ], collapsem)
}
net <- cna (cmapread, cuttof.cij=0.5,vnames=c(135:716))
save(net,file=paste0(outapp,'net_max_mod.ob'))
net = mod.select(net,"max")
save(net,file=paste0(outapp,'net_opt_mod.ob'))
netmean = mod.select(net,"mean")
save(netmean,file=paste0(outapp,'net_opt_mod_MEAN.ob'))

####Fornonlocal caca, delete above
# Rscript netob_creator_jbc.r ../../source/enswt/ultrasmall.dcd ../../source/enswt/renumber_raw.pdb ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/cijconsensus_resmat ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensus
# Rscript netob_creator_jbc.r ../../source/ensy321a/ultrasmall.dcd ../../source/ensy321a/renumber_raw.pdb ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/cijconsensus_resmat ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensus

# #run this thereafter

# python reprsentingnetob_31mar2020.py ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensusnet_max_mod.ob cij ./ ./ 1
# python reprsentingnetob_31mar2020.py ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensusnet_opt_mod.ob cij ./ ./ 1
# python reprsentingnetob_31mar2020.py ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensusnet_opt_mod_MEAN.ob cij ./ ./ 1
# python reprsentingnetob_31mar2020.py ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/amendmentconsensusnet_max_mod.ob cij ./ ./ 0
# python reprsentingnetob_31mar2020.py ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/amendmentconsensusnet_opt_mod.ob cij ./ ./ 0
# python reprsentingnetob_31mar2020.py ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/amendmentconsensusnet_opt_mod_MEAN.ob cij ./ ./ 0


# python reprsentingnetob_31mar2020.py ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensusnet_max_mod.ob cij ./ ./ 1
# python reprsentingnetob_31mar2020.py ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensusnet_opt_mod.ob cij ./ ./ 1
# python reprsentingnetob_31mar2020.py ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensusnet_opt_mod_MEAN.ob cij ./ ./ 1
# python reprsentingnetob_31mar2020.py ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/amendmentconsensusnet_max_mod.ob cij ./ ./ 0
# python reprsentingnetob_31mar2020.py ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/amendmentconsensusnet_opt_mod.ob cij ./ ./ 0
# python reprsentingnetob_31mar2020.py ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/amendmentconsensusnet_opt_mod_MEAN.ob cij ./ ./ 0

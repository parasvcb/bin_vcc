args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Please supply follwing arguements 1 cijob 2 cmapob 3 cuttoff 4 output_append", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  cijob= args[1]
  cmapob = args[2]
  cutt=as.double(args[3])
  outapp=args[4]
}
library(bio3d)
library(igraph)
print(cmapob)
load(cmapob)
print ("cmaploaded")
load(cijob)

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
net <- cna (cij, cm = cmapread, cuttof.cij=0.6,vnames=c(135:716))
save(net,file=paste0(outapp,'net_max_mod.ob'))
net = mod.select(net,"max")
save(net,file=paste0(outapp,'net_opt_mod.ob'))
netmean = mod.select(net,"mean")
save(netmean,file=paste0(outapp,'net_opt_mod_MEAN.ob'))
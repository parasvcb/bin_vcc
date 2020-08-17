
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Please supply arguements dcd and pdbfile and output_append  contactfile ", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  dcdsmall=args[1]
  pdbf=args[2]
  outapp=args[3]
  contactfile=args[4]
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
cmapread=read.table(file=contactfile,sep=",")
cmapread=as.matrix(sapply(cmapread, as.numeric))
rownames(cmapread)<-cij[,0]
colnames(cmapread)<-cij[0,]

net <- cna (cmapread, cuttof.cij=1,vnames=c(135:716),minus.log=FALSE)
mod.select <- function(x, thres=0.1) {
   remodel <- community.tree(x, rescale = TRUE)
   n.max = length(unique(x$communities$membership))
   ind.max = which(remodel$num.of.comms == n.max)
   v = remodel$modularity[length(remodel$modularity):ind.max]
   v = rev(diff(v))
   fa = which(v>=thres)[1] - 1
   ncomm = ifelse(is.na(fa), min(remodel$num.of.comms), n.max - fa)
   print (ncomm)
   ind <- which(remodel$num.of.comms == ncomm)
   network.amendment(x, remodel$tree[ind, ])
}
save(net,file=paste0(outapp,'net_max_mod.ob'))
net = mod.select(net)
save(net,file=paste0(outapp,'net_opt_mod.ob'))

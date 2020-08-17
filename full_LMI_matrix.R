args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Please supply two arguements dcd and pdbfile and output_append and numbe rof cores", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  dcdf=args[1]
  pdbf=args[2]
  outapp=args[3]
  cores=as.numeric(args[4])
}


library(bio3d)
library(igraph)
dcd <- read.dcd(dcdf)
pdb <- read.pdb(pdbf)

ca.inds <- atom.select(pdb, string="protein")

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
print ("Here")
cij<-dccm(xyz[,ca.inds$xyz], method="pearson",ncore=cores)
print ("out")
write.table(cij, file=paste0(outapp,"mymatrixcij_pearson.txt"),  sep = '\t')
save (cij,file=paste0(outapp,"CIJ_pearson.ob"))
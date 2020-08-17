args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Please supply two arguements dcd and pdbfile and output_append", call.=FALSE)
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

ca.inds <- atom.select(pdb, elety="CA")

xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
save (xyz,file=paste0(outapp,"XYZob_CA.ob"))

rf <- rmsf(xyz[,ca.inds$xyz])
rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
rdfrompdb=rmsd(pdb$xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
#df = data.frame(rmsf=rf,resid=c(135:716),grp='whole_protein',protgrp='_all')
df = data.frame(rmsf=rf,resid=c(136:716,136:716),grp='whole_protein',protgrp='_all')
dfr = data.frame(rmsd=rd,frames=c(1: dim(dcd)[1]),grp='whole_protein',protgrp='_all')
dfrfrompdb = data.frame(rmsd=rdfrompdb,frames=c(1: dim(dcd)[1]),grp='whole_protein',protgrp='_all')
write.table(df,paste0(outapp,"_rmsf.txt") ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
write.table(dfr,paste0(outapp,"_rmsd_firstframe.txt") ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
write.table(dfrfrompdb,paste0(outapp,"_rmsd_frompdb.txt") ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)

print("xyzdone")
# cij<-dccm(xyz[,ca.inds$xyz])
# print("cijdone")
# write.table(cij, file=paste0(outapp,"mymatrixcij_CA.txt"),  sep = '\t')
# print ("is")
# pdf(paste0(outapp,'dccm_CA.pdf'))
# plot(cij,sse=stride(pdb,exefile='stride',resno=TRUE))
# #plot(cij)
# save (cij,file=paste0(outapp,"CIJ_CA.ob"))

# dev.off()

pdf(paste0(outapp,'pca_res_CA.pdf'))
pc <-pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=bwr.colors(nrow(xyz)) )
dev.off()

print (pc)
save (pc,file=paste0(outapp,"pc_CA.ob"))
p1 <-mktrj.pca(pc, pc=1, b=pc$au[,1], file=paste0(outapp,"pc1_CA.pdb"))
p2 <-mktrj.pca(pc, pc=2,b=pc$au[,2], file=paste0(outapp,"pc2_CA.pdb"))
p3 <-mktrj.pca(pc, pc=3, b=pc$au[,3], file=paste0(outapp,"pc3_CA.pdb"))
p4 <-mktrj.pca(pc, pc=4, b=pc$au[,4], file=paste0(outapp,"pc4_CA.pdb"))
p5 <-mktrj.pca(pc, pc=5,b=pc$au[,5], file=paste0(outapp,"pc5_CA.pdb"))
p6 <-mktrj.pca(pc, pc=6, b=pc$au[,6], file=paste0(outapp,"pc6_CA.pdb"))

print ('pcaconformations retrieved, time for residue wise contributions to pca')
pdf(paste0(outapp,'pca_residue_contribution.pdf'))
plot.bio3d(pc$au[,1],ylab='PC1,2,3',xlab='residuepos',type='l')
points(pc$au[,2],type='l',col='blue')
points(pc$au[,3],type='l',col='green')
dev.off()
write.table(pc$au[,1],paste0(outapp,"PC1resconCA.txt"))
write.table(pc$au[,2],paste0(outapp,"PC2resconCA.txt"))
write.table(pc$au[,3],paste0(outapp,"PC3resconCA.txt"))
write.table(pc$au[,4],paste0(outapp,"PC4resconCA.txt"))
write.table(pc$au[,5],paste0(outapp,"PC5resconCA.txt"))
write.table(pc$au[,6],paste0(outapp,"PC6resconCA.txt"))

args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Please supply follwing arguements 1 ultrsmalldcd 2 cijob 3 output_append 4 contactfile", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  dcdsmall=args[1]
  cijob=args[2]
outapp=args[3]
contactfile=args[4]
}
print ('^^^^^^^contactfile')
print (args[4])
library(bio3d)
library(igraph)

cmapread=read.table(file=contactfile,sep=",")
cmapread=as.matrix(sapply(cmapread, as.numeric))
#it should be of similar format to cmapmat.txt and in matrix form, converted from df

dummydcd<-read.dcd(dcdsmall)
#2frame dcd ob

cmapdummy <- cmap(dummydcd, dcut = 4.5, scut = 0, pcut = 0.75, mask.lower = FALSE)
#creating dummy cmap object for faster processing

rownames(cmapread)<-cmapdummy[,0]

colnames(cmapread)<-cmapdummy[0,]

save(cmapread,file=paste0(outapp,"_CMAP.ob"))

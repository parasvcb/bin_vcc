library(bio3d)
library(igraph)
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("Please supply two arguements netob, and output_append and pdb and cijob", call.=FALSE)
} else if (length(args)==4) {
  # default output file
  netob=args[1]
  outapp=args[2]
  pdbf=args[3]
  obcij=args[4]
}
pdb<-read.pdb(pdbf)
func_netob <- function(netob ,outapp1) {
	summary_netob <- summary (netob)
	cat(capture.output(print(summary_netob), file=paste0(outapp1,"communities_summary.txt")))
	pdf(paste0(outapp1,'nodeplot.pdf'))
	plot( netob, pdb)
	dev.off()
	print ("2")
	pdf(paste0(outapp1,'fullplot.pdf'))
	plot( netob, pdb, full = TRUE, vertex.label.cex=0.7)
	dev.off()
	pdf(paste0(outapp1,'_and_dccmcommunities.pdf'))
	plot.dccm(cij, margin.segments =  netob$communities$membership, main="",sse=stride(pdb,exefile='stride',resno=TRUE))
    #plot.dccm(cij, margin.segments =  netob$communities$membership, main="")
	dev.off()
	pdf(paste0(outapp1,'and_dccm_filtered_node_objects.pdf'))
	plot.dccm( netob$cij, margin.segments =  netob$communities$membership, main="", sse=stride(pdb,exefile='stride',resno=TRUE))
    #plot.dccm( netob$cij, margin.segments =  netob$communities$membership, main="")
	dev.off()
	node.betweenness <- betweenness( netob$network)
	write.table(node.betweenness,paste0(outapp1,"betweencent.txt"))
	vmd.cna (netob, pdb, launch=FALSE, vmdfile = paste0(outapp1,"network.vmd"), pdbfile = paste0(outapp1, "network.pdb"))
	pdf(paste0(outapp1,'betweennenssplot.pdf'))
	plot(node.betweenness, xlab="Residue No", ylab="Centrality", type="h")
	dev.off()
	#write.pdb(pdb, b=normalize.vector(node.betweenness), file=paste0(outapp1,"betweenessCentasB.pdb"))
	pdf(paste0(outapp1,'dccmincommunitiescorr.pdf'))
	plot( netob$community.cij)
	dev.off()
    cat(capture.output(print (max( netob$communities$modularity))), file=paste0(outapp1,"max_modularity.txt"))
	write.table( netob$cij, quote=FALSE, row.names=FALSE, col.names=FALSE, file=paste0(outapp1,"adjmat.txt"))
	cat(capture.output(print(V( netob$network)$color), file=paste0(outapp1,"colorcodes_communities.txt")))
}

load(netob)
load(obcij)
func_netob(net,paste0(outapp))

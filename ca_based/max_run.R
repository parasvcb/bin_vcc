
        library(bio3d)
        library(igraph)
        load('../../derived/ensy321a/general/CIJ_CA.ob')
        load('../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.7cij/consensusnet_max_mod.ob')
        node.num <- max(net$communities$membership)
        inds <- pairwise(node.num)
        for(i in 1:nrow(inds)) {
                comms.1.inds <- which(net$communities$membership==inds[i,1])
                comms.2.inds <- which(net$communities$membership==inds[i,2])
                if (length(comms.1.inds)>2 && length(comms.2.inds)>2){
                submatrix <- net$cij[comms.1.inds, comms.2.inds]
                nonzeroindex=which(submatrix!=0, arr.ind = T)
                print (paste("
com1:",inds[i,1],"com2:",inds[i,2],"members:",sum(colSums(submatrix != 0))))
                respaircount=1
                df=as.data.frame(nonzeroindex)
                #break
                mylist <- c()
                for (j in 1:nrow(df)){
                    resind1=comms.1.inds[df[j,'row']]
                    resind2=comms.2.inds[df[j,'col']]
                    #resinds are cij indexes, for residues add 134 to them 
                    cijval=cij[resind1,resind2]
                    print (paste("-->For compair: ",inds[i,1],inds[i,2]," and residpair: ",resind1+134,resind2+134,"cij value was:",cijval))
                }
            }
        }
        
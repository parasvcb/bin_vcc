
        library(bio3d)
        library(igraph)
        load('../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.7cij/consensusnet_max_mod.ob')
        df=as_data_frame(nnet$community.network,what='both')
        print (summary(nnet))
        write.table(df$edges,'maximum_modmaxedgeweight' ,sep = "	" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
        
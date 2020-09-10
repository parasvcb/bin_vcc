netamendstring='''
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
'''
import subprocess
import os
import sys
import re
if len(sys.argv) != 4:
    print("Please enter correct cmd 1:netob 2:cijmatrix 3:outappend")
    sys.exit()

prog, netob, cijob, outputappend = sys.argv

def writehas(has,out):
    keys=list(has.keys())
    keys.sort()    
    for key in keys:
        with open (out+'key','w') as fin:
            for partitioning in has[key]:
                fin.write(" -->For partitioning :%s"%(partitioning))
                fin.write("Residuedetails:%s\n"%(has[key][partitioning]['resdet']))
                fin.write("Edgesdetails:%s\n\n\n"%(has[key][partitioning]['edges']))
            fin.write("communityMembers:%s"%(has[key]['comms']))

def getoutputofresidue(text,residuefile):
    flag= True if re.search(r'id\s+size\s+members',text) else False
    commembers = text.split("\n\n")[-2].split("members")[-1]
    has_members = {}
    for i in [i for i in commembers.split("\n") if len(i) > 1]:
        ele = i.split()
        com = int(ele[0])
        memstring = " ".join(ele[1:]) if not flag else " ".join(ele[3:])
        ref1 = re.sub('c', '', memstring)
        ref2 = re.sub(r'\(', '', ref1)
        ref3 = re.sub(r'\)', '', ref2)
        members = ref3.split(",")
        templis = []
        for j in members:
            if ":" in j:
                st, end = list(map(int, j.split(":")))
                templis += list(range(st, end+1))
            else:
                templis += [int(j)]
        has_members[com] = templis
    keys = list(has_members.keys())
    keys.sort()
    with open(residuefile, "w") as fin:
        for i in keys:
            if len(has_members[i]) > 1:
                fin.write("%s_com %s_members\nresid %s\n" % (i,len(has_members[i]), " ".join(
                    map(str, [k+134 for k in has_members[i]]))))
        fin.write("\n\npythonlistformat")
        for i in keys:
            if len(has_members[i]) > 1:
                fin.write("%s_com=[%s]\n" % (i, ",".join(
                    map(str, [k+134 for k in has_members[i]]))))
    return has_members

def give_details(megahas,members,netob,cijob,special=False,netname=False):
    #bydefault netname will be net, if sepcial is true netname can be net,netamend
    stringdefault='''
        library(bio3d)
        library(igraph)
        load('%s')
        load('%s')
        node.num <- max(%s$communities$membership)
        inds <- pairwise(node.num)
        for(i in 1:nrow(inds)) {
                comms.1.inds <- which(%s$communities$membership==inds[i,1])
                comms.2.inds <- which(%s$communities$membership==inds[i,2])
                if (length(comms.1.inds)>2 && length(comms.2.inds)>2){
                submatrix <- %s$cij[comms.1.inds, comms.2.inds]
                nonzeroindex=which(submatrix!=0, arr.ind = T)
                print (paste("\ncom1:",inds[i,1],"com2:",inds[i,2],"members:",sum(colSums(submatrix != 0))))
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
        ''' % (cijob,netob,netname,netname,netname,netname)
    moredetail=['max','mean']
    if members not in megahas:
        megahas[members]={}
    for det in moredetail:
        #print (netob,members,moredetail,members,moredetail,'nnet.ob')
        stringtype1 = '''
        library(bio3d)
        library(igraph)
        load('%s')
        %s #netamend string
        tree <- community.tree(net, rescale=TRUE)
        memb.k <- tree$tree[ tree$num.of.comms == %s, ]
        nnet <- network.amendment(net, memb.k,'%s')#mean,median tag
        node.num <- max(nnet$communities$membership)
        inds <- pairwise(node.num)
        for(i in 1:nrow(inds)) {
                comms.1.inds <- which(nnet$communities$membership==inds[i,1])
                comms.2.inds <- which(nnet$communities$membership==inds[i,2])
                if (length(comms.1.inds)>2 && length(comms.2.inds)>2){
                submatrix <- nnet$cij[comms.1.inds, comms.2.inds]
                nonzeroindex=which(submatrix!=0, arr.ind = T)
                print (paste("\ncom1:",inds[i,1],"com2:",inds[i,2],"members:",sum(colSums(submatrix != 0))))
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
        save(nnet,file='%s')
        ''' % (netob,netamendstring,members,det,members+det+'nnet.ob')
        string = stringdefault if special else stringtype1
        modnet=netob if special else members+det+'nnet.ob'
        #has stored modified object from whom residue wise details will be extracyed
        with open('%s_run.R'%det, 'w') as fin:
            fin.write("%s" % string)
        rest = subprocess.check_output(["Rscript", '%s_run.R'%det]).decode('utf-8')
        string2='''
        library(bio3d)
        library(igraph)
        load('%s')
        df=as_data_frame(nnet$community.network,what='both')
        print (summary(nnet))
        write.table(df$edges,'%s' ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
        '''%(modnet,members+det+'edgeweight') 
        with open('%s.R'%det, 'w') as fin:
            fin.write("%s" % string2)
        res = subprocess.check_output(["Rscript", '%s.R'%det]).decode('utf-8')
        has=getoutputofresidue(res,members+'xez_residuecommunitynumbering')
        
        with open (members+det+'edgeweight') as fin:
            dat=fin.read()
        if det not in megahas[members]:
            megahas[members][det]={'resdet':rest,'edges':dat}
        os.remove(members+det+'edgeweight')
        os.remove('%s.R'%det)
        os.remove('%s_run.R'%det)
    with open (members+'xez_residuecommunitynumbering') as fin:
        dat=fin.read()
    megahas[members]['comms']=dat
    os.remove(members+'xez_residuecommunitynumbering')
    return megahas

def temproutscript(netob,cijob,outputappend):
    '''
    set number of max variables, if the count is 13 then fine, else design 4 variables,
    split at 13 (based on previous), 4 (if not optimal), 6
    each of them will undergo max, mean, median and trimmed mean,
    communities all, vmd resno, python lists, tables having their max mean trimmed nmean mediana and members count and number of edges
    '''
    stringtest = '''
    library(bio3d)
    library(igraph)
    load('%s')
    node.num <- max(net$communities$membership)
    print (node.num)
    '''%netob
    with open("test_remove.r", 'w') as fin:
        fin.write("%s" % stringtest)
    res = subprocess.check_output(["Rscript", "test_remove.r"]).decode('utf-8')
    maxmembers = int(res.split()[-1])
    number_communities=[4,6,13] if maxmembers!=13 else [4,6]
    hasdetails={}
    os.remove('test_remove.r')
    #netob will be maximum
    #optob below will be optimum modulaorty counterpart
    #optobmod is above and its modified structre
    optob=re.sub(r'max','opt',netob)
    optobmod="/".join(optob.split("/")[:-1])+'amendment'+optob.split("/")[-1]
    hasdetails= give_details(hasdetails,'maximum_mod', netob,cijob, special=True,netname='net')
    hasdetails= give_details(hasdetails,'optimum_mod', optob,cijob, special=True,netname='net')
    hasdetails= give_details(hasdetails,'optimum_mod_modified', optobmod,cijob, special=True,netname='netamend')
    for coms in number_communities:
        hasdetails= give_details(hasdetails,coms, netob,cijob)
    writehas(hasdetails,outputappend)

#prog, netob, cijob, outputappend,temp = sys.argv
temproutscript(netob,cijob,outputappend)

#python netobrepresentation_nextlevel.py ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.7cij/consensusnet_max_mod.ob ../../derived/ensy321a/general/CIJ_CA.ob  ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.7cij/details_

import subprocess
import os
import sys
import re
if len(sys.argv) != 6:
    print("Please enter correct cmd 1:netob 2:type(cij or dcna) 3: inputprobfileWT 4: inputprobfileMUT (location for dcna and ./ ./ for cij) 5. neednewmembershipaspermembranebindingmotif ( 1,yes, 0 no")
    sys.exit()

prog, netob, typenet, probfilewt, probfilemt, amendment = sys.argv
amendment = int(amendment)
outputappend = ".".join(netob.split(".")[:-1])
'''
1                                                                                                                             1
2                                                                                                                             2
3                                                                                                                     c(3, 161)
4                                                                                     c(4, 57:60, 67:69, 155, 157:160, 175:178)
5  c(5:12, 50:56, 70:76, 78, 80:81, 123:136, 139, 152:154, 156, 181:183, 198, 200:204, 252:259, 261:269, 309:318, 320:321, 323)
6      c(13:19, 22, 41:49, 77, 79, 82:90, 92, 115:122, 137:138, 140:151, 184:197, 205:212, 239:251, 270:274, 276, 303:308, 490)
7                                             c(20:21, 23:40, 91, 93:99, 107:114, 213:227, 229:236, 238, 275, 277:287, 290:302)
8                                                                                                                            61
'''


def temproutscript(arg, netob, filename):
    string = '''
    library(bio3d)
    library(igraph)
    load('%s')
    df=as_data_frame(%s$community.network,what='both')
    write.table(df$edges,'%s' ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
    print (summary(%s))
    ''' % (netob, arg, outputappend+".edges", arg)
    with open(filename, 'w') as fin:
        fin.write("%s" % string)
    res = subprocess.check_output(["Rscript", filename]).decode('utf-8')
    # print(res)
    flag = True if re.search(r'id\s+size\s+members', res) else False
    commembers = res.split("\n\n")[-2].split("members")[-1]

    has_members = {}
    # colorcode=[0,1,4,7,9,10,11,16]
    for i in [i for i in commembers.split("\n") if len(i) > 1]:
        #print ("i",i)
        ele = i.split()
        com = int(ele[0])
        memstring = " ".join(ele[1:]) if not flag else " ".join(ele[3:])
        ref1 = re.sub('c', '', memstring)
        ref2 = re.sub(r'\(', '', ref1)
        ref3 = re.sub(r'\)', '', ref2)
        members = ref3.split(",")
        #print (members)
        templis = []
        for j in members:
            #print ("j",j)
            if ":" in j:
                st, end = list(map(int, j.split(":")))
                templis += list(range(st, end+1))
            else:
                templis += [int(j)]
        has_members[com] = templis
    keys = list(has_members.keys())
    keys.sort()
    # print(outputappend+".vmd_resno")
    with open(outputappend+".vmd_resno", "w") as fin:
        for i in keys:
            if len(has_members[i]) > 2:
                fin.write("%s_com\nresid %s\n" % (i, " ".join(
                    map(str, [k+134 for k in has_members[i]]))))
    return has_members


def probcounter(filename):
    has_prob = {}
    with open(filename) as fin:
        dat = [i for i in fin.read().split("\n")[1:] if len(i) > 0]
    for i in dat:
        ele = i.split("\t")
        r1, r2 = map(int, ele[0][1:-1].split(","))
        pair = (r1+135, r2+135)
        prob = float(ele[1])
        has_prob[pair] = prob
    return has_prob


def contact_prob(has, pfilewt, pfilemt):
    haswt = probcounter(pfilewt)
    hasmt = probcounter(pfilemt)
    with open(outputappend+"dcna_communityedges.edges", 'w') as fout:
        fout.write(
            "g1\tg2\tdiff(wt-mut)(intsxn)\tpairs(intsxn)\tnormalize(intsxn)\tdiff(wt-mut)(union)\tpairs(union)\tnormalize(union)\n")
        communities = list(has.keys())
        for i in range(0, len(communities)):
            for j in range(i+1, len(communities)):
                com1 = communities[i]
                com2 = communities[j]
                if len(has[com1]) > 2 and len(has[com2]) > 2:
                    # if [com1, com2] == [1, 3]:
                    #print(com1, com2)
                    su_union = []
                    su_intersection = []
                    wtcountintersec = []
                    mtcountintersec = []
                    wtcountunion = []
                    mtcountunion = []
                    # print(has)
                    for k in has[com1]:
                        k = k+134
                        for m in has[com2]:
                            m = m+134
                            pair = (k, m) if k < m else (m, k)
                            # if pair in haswt and pair in hasmt and 20 <= haswt[pair] <= 90 and 20 <= hasmt[pair] <= 90:
                            #     flag = 'intersection'
                            #     wtcount = haswt[pair]
                            #     mtcount = hasmt[pair]
                            # elif pair in haswt and pair not in hasmt and 20 <= haswt[pair] <= 90:
                            #     flag = 'union'
                            #     wtcount = haswt[pair]
                            #     mtcount = 0
                            # elif pair in hasmt and pair not in haswt and 20 <= hasmt[pair] <= 90:
                            #     flag = 'union'
                            #     mtcount = hasmt[pair]
                            #     wtcount = 0
                            # else:
                            #     flag = False
                            # different approaches different papers
                            if pair in haswt and pair in hasmt and haswt[pair] > 0 and hasmt[pair] > 0:
                                flag = 'intersection'
                                wtcount = haswt[pair]
                                mtcount = hasmt[pair]
                            elif pair in haswt and pair not in hasmt and haswt[pair] > 0:
                                flag = 'union'
                                wtcount = haswt[pair]
                                mtcount = 0
                            elif pair in hasmt and pair not in haswt and hasmt[pair] > 0:
                                flag = 'union'
                                mtcount = hasmt[pair]
                                wtcount = 0
                            else:
                                flag = False
                            if flag and abs(k-m) > 2:
                                # abs condition is for atleast 2 residues apart (noncivalent interactions)
                                # print(pair, wtcount, mtcount)
                                if flag == 'intersection':
                                    wtcountintersec += [wtcount]
                                    mtcountintersec += [mtcount]
                                    su_intersection += [wtcount-mtcount]
                                su_union += [wtcount-mtcount]
                                wtcountunion += [wtcount]
                                mtcountunion += [mtcount]
                    norm_union = sum(su_union) / \
                        len(su_union) if su_union else 0
                    groupcountunion = (sum(wtcountunion)-sum(mtcountunion))
                    groupcountintersec = (
                        sum(wtcountintersec)-sum(mtcountintersec))
                    normgintersec = groupcountintersec / \
                        len(wtcountintersec) if len(wtcountintersec) > 0 else 0
                    normgunion = groupcountunion / \
                        len(wtcountunion) if len(wtcountunion) > 0 else 0
                    norm_intersection = sum(
                        su_intersection)/len(su_intersection) if su_intersection else 0
                    fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (com1, com2, sum(su_intersection)/100, len(
                        su_intersection), norm_intersection/100, sum(su_union)/100, len(su_union), norm_union/100))
                    # fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgroupaddition\n" % (com1, com2, groupcountintersec, len(
                    #    wtcountintersec), normgintersec, groupcountunion, len(wtcountunion), normgunion))


def amendlist(has, filename):
    # print(has)
    resid_vector = {}
    for key in has:
        if len(has[key]) > 1:
            for resid in has[key]:
                resid_vector[resid] = key
    needed_mem = set(list(resid_vector.values()))
    # needed mem will be a list of `communities`
    needed_com = 1
    their_hash = {}
    for j in needed_mem:
        their_hash[j] = needed_com
        needed_com += 1
    motif = needed_com
    # their_hash will have commids has keys and serial ids as as new keys startring from 1, motif got new membership
    # new membership
    needed_com += 1
    with open(filename, 'w') as fin:
        # fin.write("Membership\n")
        for j in range(1, 583):
            # need total of 582 elemnets:
            if j+134 in range(292, 311):
                fin.write("%s\n" % motif)
            else:
                if j in resid_vector:
                    fin.write("%s\n" % (their_hash[resid_vector[j]]))
                else:
                    fin.write("%s\n" % needed_com)
                    needed_com += 1
    return


def amendment_func(netob, vectorlist, arg, save_loc, meanormax):
    string = '''
    library(bio3d)
    library(igraph)
    load('%s')
    a=read.csv('%s',header=FALSE)
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

    netamend <- network.amendment(%s, a$V1, "%s")
    save(netamend,file='%s')
    ''' % (netob, vectorlist, arg, meanormax, save_loc)
    with open("amednnet.R", 'w') as fin:
        fin.write("%s" % string)
    res = subprocess.check_output(["Rscript", "amednnet.R"]).decode('utf-8')
    return


nettypestring = netob.split("/")[-1]
argobname = 'netmean' if 'MEAN' in nettypestring else 'net'
print('amendment' in netob)
print(amendment)
if 'amendment' in netob and amendment == 0:
    print("**yes")
    argobname = "netamend"

has_members = temproutscript(argobname, netob, "tempscript.R")
if typenet == "dcna":
    contact_prob(has_members, probfilewt, probfilemt)
if amendment:
    amendlist(has_members, outputappend+".amendedvectorlist")
    newfile = "/".join(netob.split("/")[:-1])+"/amendment"+netob.split("/")[-1]
    print(newfile)
    meanormax = "mean" if "MEAN" in netob else "max"
    amendment_func(netob, outputappend +
                   ".amendedvectorlist", argobname, newfile, meanormax)

import sys
import pickle
import subprocess
import re
import os
import mdtraj
import numpy as np
from collections import Counter


def callscript(text, filename, mode='R'):
    if mode == 'R':
        with open(filename, 'w') as fout:
            fout.write(text)
        res = subprocess.check_output(['Rscript', filename]).decode('utf-8')
        os.remove(filename)
        return res


#############################
# CONTACTS AREA
# PROCESSING CONTACTS AND MASTERFILE
#############################

def matrix_returner_single_cuttoff(numpy_has_file, distcutt, range_cutt, upp=False, low=False):
    with open(numpy_has_file, 'rb') as fin:
        numpy_has = pickle.load(fin)
    # first rep will have 0:10000, seoc nd will have 10000:20000, 20000:30000
    # rangecutt will be less than 1,
    # distcutt should be in angstrom
    range_cutt = range_cutt if range_cutt <= 1 else range_cutt/100
    probhas = {}
    keys = list(numpy_has.keys())
    if upp == False:
        upp = len(numpy_has[keys[0]])
    if low == False:
        low = 0
    last = max([i[1] for i in keys])
    mat = np.zeros((last+1, last+1))
    middfrac = ((upp-low)//2)+low
    for residpair in numpy_has:
        valuesfull = np.where(numpy_has[residpair][low:upp] <= distcutt)
        valuesfhalf = np.where(numpy_has[residpair][low:middfrac] <= distcutt)
        valuesshalf = np.where(numpy_has[residpair][middfrac:upp] <= distcutt)
        resultfull = len(valuesfull[0])/(upp-low)
        resultfhalf = len(valuesfhalf[0])/(middfrac-low)
        resultshalf = len(valuesshalf[0])/(upp-middfrac)
        val = 1 if resultfull >= range_cutt else 0
        mat[residpair[0], residpair[1]] = val
        mat[residpair[1], residpair[0]] = val
        probhas[residpair] = [resultfull, resultfhalf, resultshalf]
    # matrix is computed
    return mat, probhas


def contacts_computer(has, dcd, pdb, contype=1):
    '''
    residue number will occasonally sgtart from 0, add 135 to them later
    has={
        (first,second):np.array([].dtype='float16')
    }
    '''
    print("for dcd %s" % dcd)

    t = mdtraj.load_dcd(dcd, top=pdb)
    print("Loaded")
    if int(contype) == 0:
        print('computing closest heavy atom')
        con = mdtraj.compute_contacts(
            t, contacts='all', scheme='closest-heavy', ignore_nonprotein=True)
    else:
        print('computing CA-CA distance')
        con = mdtraj.compute_contacts(
            t, contacts='all', scheme='ca', ignore_nonprotein=True)
    print("computed")
    object0, object1 = con
    if len(has) == 0:
        for i in object1:
            has[tuple(i)] = np.array([], dtype='float16')
    reslis = object0[0]
    for respairind in range(0, len(reslis)):
        distance_in_frames = object0[:, respairind]
        correspondingresiduepair = tuple(object1[respairind])
        has[correspondingresiduepair] = np.append(
            has[correspondingresiduepair], np.round(distance_in_frames*10, 4))
    return has


def trajsplitter(dcd, pdb, splitfac):
    tet = subprocess.check_output(["catdcd", dcd])
    tet = tet.decode("utf-8")
    frames = int(re.search(r'Read \d+ frames', tet).group().split()[1])
    steps = frames//splitfac
    rem = frames % splitfac
    tlis = []
    for i in range(0, steps):
        addfac = 0 if i == 0 else 1
        tlis += [((i*splitfac)+addfac, (i+1)*splitfac)]
        i = i+1
        if rem:
            tlis += [((i*splitfac)+addfac, (i*splitfac)+rem)]
    for ind, val in enumerate(tlis):
        createdcd = subprocess.check_output(
            ["catdcd", "-o", "%s.tempdcd" % ind, "-first", " %s" % val[0], "-last", " %s" % val[1], dcd])

    return tlis


def horse(dcd, pdb, outname, contype=1):
    # contype=1 means calcluate caca non local distances, 0 means closet heavy CACA nonlocal distance
    tlis = trajsplitter(dcd, pdb, 1000)
    append = "contacts_ca-ca_nonlocal.hash" if int(
        contype) else "contacts_closest_heavy_nonlocal.hash"
    has_mega = {}
    # this should be offloaded and added if needed somewhere,a nd can be reuploaded
    for ind, val in enumerate(tlis):
        has_mega = contacts_computer(
            has_mega, "%s.tempdcd" % ind, pdb, contype)
        os.remove("%s.tempdcd" % ind)
    with open(os.path.join(outname, append), "wb") as fin:
        pickle.dump(has_mega, fin, protocol=pickle.HIGHEST_PROTOCOL)


#############################
# GETTING CIJ
# SMALL PROGRAMS
# CREATING CONSENSUS MATRIX
# GETTING SPARSENESS INFORMATION
# GENERATING GENERAL INFORMATION
#############################
def write_matrix(matrix, outfile):
    with open(outfile+"_resmat", "w") as fin:
        firststring = ''
        for i in range(0, len(matrix)):
            firststring += '"V%s",' % (i+1)
        firststring = firststring[:-1]
        fin.write("%s\n" % firststring)
        for ind, val in enumerate(matrix):
            linstr = ''
            for j in val:
                linstr += "%s," % j
            fin.write('"%s",%s\n' % (ind+1, linstr[:-1]))
    np.save(outfile+"_.npyob", matrix, allow_pickle=True)


def writeprobfile(probhas, outfile):
    #
    keys = list(probhas.keys())
    keys.sort()
    outfile = re.sub(r'\d+pcnt.*', '_prob', outfile)
    with open(outfile, "w") as fin:
        fin.write(
            "Residpair\tP(totalframes)\tP(first50frames)\tP(second50frames)\tsubP(first-second)\n")
        for i in keys:
            if probhas[i][0] > 0:
                fullval, fhalf, shalf = probhas[i]
                sub = (fhalf-shalf)*100
                fin.write("%s\t%s\t%s\t%s\t%s\n" %
                          (i, fullval*100, fhalf*100, shalf*100, sub))


def texttomatrix(matfile):
    with open(matfile) as fin:
        dat = fin.read()
    dat = re.sub(r'"', '', dat)
    dat = re.sub(r'V', '', dat)
    with open("testmat", 'w') as fin:
        fin.write("%s" % dat)
    a = np.loadtxt(open("testmat"), delimiter="\t", skiprows=1)
    a = [i[1:] for i in a]
    # a= [[round(j,3)] for i in a for j in i]
    os.remove('testmat')
    return np.array(a)


def callcij(dcd, pdb, outadd):
    string = '''
    library(bio3d)
    dcd <- read.dcd('%s')
    pdb <- read.pdb('%s')
    ca.inds <- atom.select(pdb, elety="CA")
    xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
    cij<-dccm(xyz[,ca.inds$xyz])
    print("cijdone")
    write.table(cij, file=file.path("%s","mymatrixcij_CA.txt"),  sep = '\t')
    print ("is")
    pdf(file.path('%s','dccm_CA.pdf'))
    plot(cij,sse=stride(pdb,exefile='stride',resno=TRUE))
    plot(cij)
    save (cij,file=file.path('%s',"CIJ_CA.ob"))
    dev.off()
    ''' % (dcd, pdb, outadd, outadd, outadd)
    callscript(string, 'tempscript.R')


def consensus_matrix(mat1, mat2, mat3, con1, con2, con3, cuttcij, outfile):
    # mat1, mat2, mat3, con1, con2, con3, cuttcij, consensusFile
    def tempreturn(val, denom):
        if val:
            return round(val/denom, 3)
        else:
            return 0
    totalpairs = len(mat1)*len(mat1)
    string = 'Total:%s\n' % (totalpairs)
    combinedmatrix = (abs(mat1)+abs(mat2)+abs(mat3))/3
    temp1 = len(np.where(combinedmatrix >= 0.1)[0])
    temp2 = len(np.where(combinedmatrix >= 0.2)[0])
    temp3 = len(np.where(combinedmatrix >= 0.3)[0])
    temp4 = len(np.where(combinedmatrix >= 0.4)[0])
    temp5 = len(np.where(combinedmatrix >= 0.5)[0])
    temp6 = len(np.where(combinedmatrix >= 0.6)[0])
    temp7 = len(np.where(combinedmatrix >= 0.7)[0])
    temp8 = len(np.where(combinedmatrix >= 0.8)[0])
    temp9 = len(np.where(combinedmatrix >= 0.9)[0])
    string += '''Total greater than 0.1 were: %s (%s)
    Total greater than 0.2 were: %s (%s)
    Total greater than 0.3 were: %s (%s)
    Total greater than 0.4 were: %s (%s)
    Total greater than 0.5 were: %s (%s)
    Total greater than 0.6 were: %s (%s)
    Total greater than 0.7 were: %s (%s)
    Total greater than 0.8 were: %s (%s)
    Total greater than 0.9 were: %s (%s)
    ''' % (temp1, tempreturn(temp1, totalpairs), temp2, tempreturn(temp2, totalpairs),
           temp3, tempreturn(temp3, totalpairs), temp4, tempreturn(
               temp4, totalpairs),
           temp5, tempreturn(temp5, totalpairs), temp6, tempreturn(
               temp6, totalpairs),
           temp7, tempreturn(temp7, totalpairs), temp8, tempreturn(
               temp8, totalpairs),
           temp9, tempreturn(temp9, totalpairs))
    string += '\nGeneral matrix sparsity analysed\n'
    hascount = {}
    hasvalues = {}
    consensuscount = 0
    distancequalified = 0
    oppsign = 0
    oppsignqualified = 0

    cons_np_cij = np.zeros((len(mat1), len(mat1)))
    minhas = {0: con1, 1: con2, 2: con3}
    for i in range(0, len(mat1)):
        for j in range(0, len(mat1[i])):
            val1 = mat1[i][j]
            val2 = mat2[i][j]
            val3 = mat3[i][j]
            tempkey = (i+135, j+135)
            cond1 = abs(val1) >= cuttcij and abs(
                val2) >= cuttcij and abs(val3) >= cuttcij
            cond2 = np.sign(val1) == np.sign(val2) == np.sign(val3)
            if cond1 and cond2:
                # avg=(val1+val2+val3)/3
                avg = (abs(val1)+abs(val2)+abs(val3))/3
                cons_np_cij[i][j] = avg
                consensuscount += 1
                hascount[tempkey] = 'consensus'
                hasvalues[tempkey] = round(avg, 3)
            else:
                if cond1 and not cond2:
                    hascount[tempkey] = 'opp sign_'
                    oppsign += 1
                # chcek if any two have greater tthan cuttcij then if its maintained in 75% or so frames
                flag = False
                if abs(val1) >= cuttcij:
                    if minhas[0][i][j] == 1:
                        cons_np_cij[i][j] = val1
                        flag = True
                if flag == False and abs(val2) >= cuttcij:
                    if minhas[1][i][j] == 1:
                        cons_np_cij[i][j] = val2
                        flag = True
                if flag == False and abs(val3) >= cuttcij:
                    if minhas[2][i][j] == 1:
                        cons_np_cij[i][j] = val3
                        flag = True
                if flag == False:
                    cons_np_cij[i][j] = 0
                hasvalues[tempkey] = round(cons_np_cij[i][j], 3) if flag else 0
                tempvalue = 'distanceSatisfied' if flag else 'unqualified'
                if tempkey not in hascount:
                    hascount[tempkey] = ''
                hascount[tempkey] += tempvalue
                if flag:
                    distancequalified += 1
                if cond1 and not cond2 and flag:
                    oppsignqualified += 1
    write_matrix(cons_np_cij, outfile)
    string += '''
    consensusPairs=%s
    oppositeSign=%s
    distanceQualified=%s
    oppSign and qualified=%s
    ''' % (consensuscount, oppsign, distancequalified, oppsignqualified)
    keys = list(hascount.keys())
    keys.sort()
    # print('is>>')
    with open(outfile+'_qualifierType.hashtext', 'w') as fout:
        fout.write("res1\tres2\tqualifier\n")
        for i in keys:
            fout.write("%s\t%s\t%s\n" % (i[0], i[1], hascount[i]))
    with open(outfile+'_valuesCij.hashtext', 'w') as fout:
        fout.write("res1\tres2\tvalue\n")
        for i in keys:
            fout.write("%s\t%s\t%s\n" % (i[0], i[1], hasvalues[i]))
    with open(outfile+'_generalInfo.text', 'w') as fout:
        fout.write("%s" % string)
    return

#############################
# CREATINGNETOB
#############################


def createnetob(ultrasmalldcd, pdbf, consensuscij, outloc):
    stringprogram = '''
    library(bio3d)
    library(igraph)
    dummydcd<-read.dcd('%s')
    pdb <- read.pdb('%s')
    ca.inds <- atom.select(pdb, elety="CA")
    xyz <- fit.xyz(fixed=pdb$xyz, mobile=dummydcd,fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
    print("xyzdone")
    cij<-dccm(xyz[,ca.inds$xyz])
    # from dummy dcd
    cmapread=read.table(file='%s',sep=",")
    cmapread=as.matrix(sapply(cmapread, as.numeric))
    rownames(cmapread)<-cij[,0]
    colnames(cmapread)<-cij[0,]

    net <- cna (cmapread, cuttof.cij=0.6,vnames=c(135:716))
    save(net,file=file.path('%s','net_max_mod.ob'))

    # Fornonlocal caca, delete above
    # Rscript netob_creator_jbc.r ../../source/enswt/ultrasmall.dcd ../../source/enswt/renumber_raw.pdb ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/cijconsensus_resmat ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensus
    # python reprsentingnetob_31mar2020.py ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/consensusnet_max_mod.ob cij ./ ./ 1
    ''' % (ultrasmalldcd, pdbf, consensuscij, outloc)
    callscript(stringprogram, 'netob.R', mode='R')

#############################
# MICROINFORMATION FROM NETOB
#############################

# GETCOMMUNITYFILEANDHASH (COMMUNITY FILE, GENERAL INFO (DF EDGES))


def get_members(netob, outapp):
    string = '''
    library(bio3d)
    library(igraph)
    load('%s')
    df=as_data_frame(%s$community.network,what='both')
    write.table(df$edges,'%s' ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
    print (summary(%s))
    ''' % (netob, 'net', outapp+".edges", 'net')
    res = callscript(string, 'temp.R')
    # print(res)
    flag = True if re.search(r'id\s+size\s+members', res) else False
    commembers = res.split("\n\n")[-2].split("members")[-1]

    has_members = {}
    # colorcode=[0,1,4,7,9,10,11,16]
    for i in [i for i in commembers.split("\n") if len(i) > 1]:
        # print ("i",i)
        ele = i.split()
        com = int(ele[0])
        memstring = " ".join(ele[1:]) if not flag else " ".join(ele[3:])
        ref1 = re.sub('c', '', memstring)
        ref2 = re.sub(r'\(', '', ref1)
        ref3 = re.sub(r'\)', '', ref2)
        members = ref3.split(",")
        # print (members)
        templis = []
        for j in members:
            # print ("j",j)
            if ":" in j:
                st, end = list(map(int, j.split(":")))
                templis += list(range(st, end+1))
            else:
                templis += [int(j)]
        has_members[com] = templis
    # print(has_members)
    keys = list(has_members.keys())
    keys.sort()
    with open(outapp+".vmd_resno", "w") as fin:
        for i in keys:
            if len(has_members[i]):
                fin.write("%s_com %s_members\nresid %s\n" % (i, len(has_members[i]), " ".join(
                    map(str, [k+134 for k in has_members[i]]))))
        fin.write("\n\npythonlistformat")
        for i in keys:
            if len(has_members[i]):
                fin.write("%s_com=[%s]\n" % (i, ",".join(
                    map(str, [k+134 for k in has_members[i]]))))
    return has_members

# PARSETHROUGH R (GENERAL)


def quickRparse(netob, cij_rep_any, cmapread, outfile):
    stringtype1 = '''
    library(bio3d)
    library(igraph)
    load('%s')
    load('%s')
    cmapread=read.table(file='%s',sep=",")
    cmapread=as.matrix(sapply(cmapread, as.numeric))
    rownames(cmapread)<-cij[,0]
    colnames(cmapread)<-cij[0,]
    cij=cmapread
    node.num <- max(net$communities$membership)
    inds <- pairwise(node.num)
    for(i in 1:nrow(inds)) {
            comms.1.inds <- which(net$communities$membership==inds[i,1])
            comms.2.inds <- which(net$communities$membership==inds[i,2])
            if (length(comms.1.inds)>1 && length(comms.2.inds)>1){
                submatrix <- net$cij[comms.1.inds, comms.2.inds]
                if (sum(colSums(submatrix != 0))>0) {
                    nonzeroindex=which(submatrix!=0, arr.ind = T)
                    print (paste("PARENT","-->ForCompair:",inds[i,1],inds[i,2],"members:",sum(colSums(submatrix != 0))),sep=":\\t")
                    respaircount=1
                    df=as.data.frame(nonzeroindex)
                    #break
                    mylist <- c()
                    for (j in 1:nrow(df)){
                        resind1=comms.1.inds[df[j,'row']]
                        resind2=comms.2.inds[df[j,'col']]
                        #resinds are cij indexes, for residues add 134 to them 
                        cijval=cij[resind1,resind2]
                        print (paste("CHILD","-->ForCompair: ",inds[i,1],inds[i,2]," and residpair: ",resind1+134,resind2+134,"cij value was:",cijval),sep=":\\t")
                    }
            }
        }
        else{
            print (paste("EMPTY","Compair:",inds[i,1],inds[i,2]),sep=":\\t")
        }

    }
    ''' % (netob, cij_rep_any, cmapread)
    res = callscript(stringtype1, 'temp.R')
    with open(outfile, 'w') as fout:
        fout.write("%s" % (res))
    return res


def immediateFetch(key, hascom, cijhas):
    members1 = hascom[int(key[0])]
    members2 = hascom[int(key[1])]
    value = []
    for i in members1:
        for j in members2:
            key = [i, j]
            key.sort()
            key = tuple(key)
            if key in cijhas:
                value += [(key[0], key[1], cijhas(key))]
    return value


def processfinaldata(deta, cijfile, hasmembership, qualifierFile, finalFile):
    # print(cijfile)
    deta = re.sub(r'"', '', deta)
    deta = [i for i in deta.split("\n")if len(i) > 0]
    # print(deta[:10])
    with open(cijfile) as fin:
        dat = [i for i in fin.read().split("\n")[1:] if len(i) > 0]
    hascij = {}
    hascommunities = {}
    for i in dat:
        r1, r2, val = i.split()
        r1 = int(r1)
        r2 = int(r2)
        val = float(val)
        key = [r1, r2]
        key.sort()
        hascij[tuple(key)] = round(val, 3)
    tempdata = []
    key = []
    for record in deta:
        record = re.sub(r'^\[\d+\]\s+', '', record)
        # print(record)
        if len(record) > 0:
            # print(record)
            ele = record.split()
            if ele[0] == 'PARENT':
                if tempdata and key:
                    hascommunities[key] = tempdata
                    tempdata = []
                key = (int(ele[2]), int(ele[3]))
            if ele[0] == 'CHILD':
                tempdata += [(ele[6], ele[7], ele[11])]
            if ele[0] == 'EMPTY':
                if tempdata and key:
                    hascommunities[key] = tempdata
                    tempdata = []
                key = (int(ele[2]), int(ele[3]))
                tempdata = immediateFetch(key, hasmembership, hascij)
                #print(record, tempdata)
    if tempdata and key:
        hascommunities[key] = tempdata
        tempdata = []
    keys = list(hascommunities.keys())
    keys.sort()
    hasqualifier = {}
    with open(qualifierFile) as fin:
        dat = [i for i in fin.read().split('\n')[1:] if len(i) > 1]
        for i in dat:
            ele = i.split()
            if ele[2] != 'unqualified':
                key = (ele[0], ele[1])
                hasqualifier[key] = ele[2]
    with open(finalFile, 'w') as fout:
        strings = []
        for i in keys:
            vals = [[abs(float(k[2])), (k[0], k[1])]
                    for k in hascommunities[i]]
            allqualifiers = [hasqualifier[k[1]] for k in vals]
            allvals = [k[0] for k in vals]
            vals.sort()
            tempstr = '''Com pair: %s and %s,TotalInterComEdges:%s, totalCorrelationVal:%s, 
            meanval:%s, maxval:%s (residpair:%s, qualifier:%s)
            For all residues %s were \n''' % (i[0], i[1], len(allvals), sum(allvals), np.mean(allvals), vals[-1][0], vals[-1][1], hasqualifier[vals[-1][1]], Counter(allqualifiers))
            fout.write(tempstr)
            strings += [tempstr]
            for j in hascommunities[i]:
                fout.write("\tresiduePair: %s and %s\tcijval:%s\n" %
                           (int(j[0]), int(j[1]), j[2]))
            fout.write("\n")
        stringsData = "\n\n".join(strings)
        fout.write('%s' % (stringsData))
    return
# PARSE R OUTPUT WITH PYTHON TO GET FINAL FILE

# get GENERAL analysis


def generalAnalysis(dcd, pdb, outapp):
    string = '''
    library(bio3d)
    library(igraph)
    dcd <- read.dcd('%s')
    pdb <- read.pdb('%s')
    ca.inds <- atom.select(pdb, elety="CA")
    xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
    rf <- rmsf(xyz[,ca.inds$xyz])
    rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
    rdfrompdb=rmsd(pdb$xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
    df = data.frame(rmsf=rf,resid=c(135:716),grp='whole_protein',protgrp='_all')
    dfr = data.frame(rmsd=rd,frames=c(1: dim(dcd)[1]),grp='whole_protein',protgrp='_all')
    dfrfrompdb = data.frame(rmsd=rdfrompdb,frames=c(1: dim(dcd)[1]),grp='whole_protein',protgrp='_all')
    write.table(df,paste0('%s',"_rmsf.txt") ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
    write.table(dfr,paste0('%s',"_rmsd_firstframe.txt") ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
    write.table(dfrfrompdb,paste0('%s',"_rmsd_frompdb.txt") ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
    ''' % (dcd, pdb, outapp, outapp, outapp)
    res = callscript(string, 'general.R')
    return

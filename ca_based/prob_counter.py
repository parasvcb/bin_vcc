# will fetch the matrices from numpy objects
# or will give the proablity ocuureence list

import numpy as np
import pickle
import sys,re
if len(sys.argv) != 7:
    print("Please enter correct arguments: 1 numpyhash 2 distancecuttoff(in ang)\
    3. rangecuttoff (singlevalue for greaterthna cutt, eg 0.75 or 0.2:0.7 for range cutt)\
    4. frames cuttoff in range like  0:10000, 10000:20000, 20000:30000 else False \
    5. outfile 6, ifprobfile_needstobewritten (1 for yes or 0 for no)")
    sys.exit()
pythonprog, numpyhasfile, distcutt, rangecutt, framescutt, outfile, probflag = sys.argv
# get the arguments

def matrix_returner_single_cuttoff(numpy_has, distcutt, range_cutt, upp=False, low=False):
    # first rep will have 0:10000, seoc nd will have 10000:20000, 20000:30000
    # rangecutt will be less than 1,
    # distcutt should be in angstrom
    range_cutt = range_cutt if range_cutt <= 1 else range_cutt/100
    probhas={}
    keys = list(numpy_has.keys())
    if upp == False:
        upp = len(numpy_has[keys[0]])
    if low == False:
        low = 0
    last = max([i[1] for i in keys])
    mat = np.zeros((last+1, last+1))
    middfrac=((upp-low)//2)+low
    count=0
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
        probhas[residpair]=[resultfull,resultfhalf,resultshalf]
    #matrix is computed
    return mat,probhas


def matrix_returner_range_cuttoff(numpy_has, distcutt, range_cuttlow, range_cuttupp, upp=False, low=False):
    
    # first rep will have 0:10000, seoc nd will have 10000:20000, 20000:30000
    range_cuttlow = range_cuttlow if range_cuttlow <= 1 else range_cuttlow/100
    range_cuttupp = range_cuttupp if range_cuttupp <= 1 else range_cuttupp/100
    probhas={}
    keys = list(numpy_has.keys())
    if upp == False:
        upp = len(numpy_has[keys[0]])
    if low == False:
        low = 0
    middfrac=((upp-low)//2)+low
    last = max([i[1] for i in keys])
    mat = np.zeros((last+1, last+1))
    for residpair in numpy_has:
        valuesfull = np.where(numpy_has[residpair][low:upp] <= distcutt)
        valuesfhalf = np.where(numpy_has[residpair][low:middfrac] <= distcutt)
        valuesshalf = np.where(numpy_has[residpair][middfrac:upp] <= distcutt)
        # (array([0, 1, 2, 3], dtype=int32),),, storedin vaues
        resultfull = len(valuesfull[0])/(upp-low)
        resultfhalf = len(valuesfhalf[0])/(middfrac-low)
        resultshalf = len(valuesshalf[0])/(upp-middfrac)
        val = 1 if range_cuttlow <= resultfull <= range_cuttupp else 0
        mat[residpair[0], residpair[1]] = val
        mat[residpair[1], residpair[0]] = val
        probhas[residpair]=[resultfull,resultfhalf,resultshalf]
    return mat,probhas


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

def writeprobfile(probhas,outfile):
    #
    keys=list(probhas.keys())
    keys.sort()
    outfile=re.sub(r'\d+pcnt.*','_prob',outfile)
    with open (outfile,"w") as fin:
        fin.write("Residpair\tP(totalframes)\tP(first50frames)\tP(second50frames)\tsubP(first-second)\n")
        for i in keys:
            if probhas[i][0]>0:
                fullval,fhalf,shalf=probhas[i]
                sub=(fhalf-shalf)*100
                fin.write("%s\t%s\t%s\t%s\t%s\n"%(i,fullval*100,fhalf*100,shalf*100,sub))
    
# pythonprog,numpyhasfile,distcutt,rangecutt,framescutt,outfile=sys.argv
distcutt = float(distcutt)
# this variable value shlud ibe in angstrom
uprange, lowrange = [False, False]
flow, fupp = [False, False] if framescutt == "False" else list(
    map(int, framescutt.split(':')))
# above value needed if partciular ensemble matrix needs to be generated, default hwole ensemble
with open(numpyhasfile, "rb") as fin:
    numpyhas = pickle.load(fin)
range_cutt_val = float(rangecutt) if ":" not in rangecutt else False
if range_cutt_val == False:
    lowrange, uprange = list(map(float, rangecutt.split(":")))
    mat,probhas = matrix_returner_range_cuttoff(
        numpyhas, distcutt, lowrange, uprange, fupp, flow)
else:
    mat,probhas = matrix_returner_single_cuttoff(
        numpyhas, distcutt, range_cutt_val, fupp, flow)

write_matrix(mat, outfile)
if int(probflag):
    writeprobfile(probhas,outfile)

# python prob_counter.py ../../source/enswt/contacts_ca-ca_nonlocal.hash 10 0.75 0:10000 ../../derived/wt_r1/general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes 1
# python prob_counter.py ../../source/enswt/contacts_ca-ca_nonlocal.hash 10 0.75 10000:20000 ../../derived/wt_r2/general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes 1
# python prob_counter.py ../../source/enswt/contacts_ca-ca_nonlocal.hash 10 0.75 20000:30000 ../../derived/wt_r3/general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes 1

# python prob_counter.py ../../source/ensy321a/contacts_ca-ca_nonlocal.hash 10 0.75 0:10000 ../../derived/y321a_r4/general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes 1
# python prob_counter.py ../../source/ensy321a/contacts_ca-ca_nonlocal.hash 10 0.75 10000:20000 ../../derived/y321a_r2/general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes 1
# python prob_counter.py ../../source/ensy321a/contacts_ca-ca_nonlocal.hash 10 0.75 20000:30000 ../../derived/y321a_r3/general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes 1



##fordcnanetworks
#  python prob_counter.py ../../source/enswt/allcontacts_closest_heavy.hash 4.5 0.90 False ../../derived/enswt/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes 1
#  python prob_counter.py ../../source/ensy321a/allcontacts_closest_heavy.hash 4.5 0.90 False ../../derived/ensy321a/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes 1
 

#  python prob_counter.py ../../source/enswt/allcontacts_closest_heavy.hash 4.5 0.90 0:10000 ../../derived/wt_r1/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes 1
#  python prob_counter.py ../../source/enswt/allcontacts_closest_heavy.hash 4.5 0.90 10000:20000 ../../derived/wt_r2/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes 1
#  python prob_counter.py ../../source/enswt/allcontacts_closest_heavy.hash 4.5 0.90 20000:30000 ../../derived/wt_r3/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes 1
 
#  python prob_counter.py ../../source/ensy321a/allcontacts_closest_heavy.hash 4.5 0.90 0:10000 ../../derived/y321a_r4/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes 1
#  python prob_counter.py ../../source/ensy321a/allcontacts_closest_heavy.hash 4.5 0.90 10000:20000 ../../derived/y321a_r2/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes 1
#  python prob_counter.py ../../source/ensy321a/allcontacts_closest_heavy.hash 4.5 0.90 20000:30000 ../../derived/y321a_r3/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes 1

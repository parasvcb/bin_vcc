# will fetch the matrices from numpy objects
# or will give the proablity ocuureence list

import numpy as np
import pickle
import sys
if len(sys.argv) != 6:
    print("Please enter correct arguments: 1 numpyhash 2 distancecuttoff(in ang)\
    3. rangecuttoff (singlevalue for greaterthna cutt, eg 0.75 or 0.2:0.7 for range cutt)\
    4. frames cuttoff in range like  0:10000, 10000:20000, 20000:30000 else False \
    5. outfile")
    sys.exit()
pythonprog, numpyhasfile, distcutt, rangecutt, framescutt, outfile = sys.argv
# get the arguments

def matrix_returner_single_cuttoff(numpy_has, distcutt, range_cutt, upp=False, low=False):
    # first rep will have 0:10000, seoc nd will have 10000:20000, 20000:30000
    # rangecutt will be less than 1,
    # distcutt should be in angstrom
    range_cutt = range_cutt if range_cutt <= 1 else range_cutt/100
    keys = list(numpy_has.keys())
    if upp == False:
        upp = len(numpy_has[keys[0]])
    if low == False:
        low = 0
    last = max([i[1] for i in keys])
    mat = np.zeros((last+1, last+1))
    count=0
    for residpair in numpy_has:
        values = np.where(numpy_has[residpair][low:upp] <= distcutt)
        result = len(values[0])/(upp-low)
        if result<=0.9 and result>=0.2:
            count+=1
        val = 1 if result >= range_cutt else 0
        mat[residpair[0], residpair[1]] = val
        mat[residpair[1], residpair[0]] = val
    #matrix is computed
    print (count,"distcutt",distcutt)
    return mat


def matrix_returner_range_cuttoff(numpy_has, distcutt, range_cuttlow, range_cuttupp, upp=False, low=False):
    
    # first rep will have 0:10000, seoc nd will have 10000:20000, 20000:30000
    range_cuttlow = range_cuttlow if range_cuttlow <= 1 else range_cuttlow/100
    range_cuttupp = range_cuttupp if range_cuttupp <= 1 else range_cuttupp/100
    
    keys = list(numpy_has.keys())
    if upp == False:
        upp = len(numpy_has[keys[0]])
    if low == False:
        low = 0
    last = max([i[1] for i in keys])
    mat = np.zeros((last+1, last+1))
    for residpair in numpy_has:
        values = np.where(numpy_has[residpair][low:upp] <= distcutt)
        # (array([0, 1, 2, 3], dtype=int32),),, storedin vaues
        result = len(values[0])/(upp-low)
        
        val = 1 if range_cuttlow <= result <= range_cuttupp else 0
        mat[residpair[0], residpair[1]] = val
        mat[residpair[1], residpair[0]] = val
    return mat


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
    mat = matrix_returner_range_cuttoff(
        numpyhas, distcutt, lowrange, uprange, fupp, flow)
else:
    mat = matrix_returner_single_cuttoff(
        numpyhas, distcutt, range_cutt_val, fupp, flow)

write_matrix(mat, outfile)

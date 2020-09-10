import numpy as np
import sys
import re

if len(sys.argv) != 7:
    print("Please enter correct cmd rep1add rep2add rep3add cijcutt postappendname(CA or fullatom) consensusmatrixname")
    sys.exit()
prog, add1, add2, add3, cuttcij, postappend, resfile = sys.argv
cuttcij = float(cuttcij)


def texttomatrix(matfile):
    with open(matfile) as fin:
        dat = fin.read()
    dat = re.sub(r'"', '', dat)
    dat = re.sub(r'V', '', dat)
    with open("testmat", 'w') as fin:
        fin.write("%s" % dat)
    a = np.loadtxt(open("testmat"), delimiter="\t", skiprows=1)
    a = [i[1:] for i in a]
    #a= [[round(j,3)] for i in a for j in i]
    return np.array(a)


mat1contact = add1 + \
    "general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes_.npyob.npy"
mat2contact = add2 + \
    "general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes_.npyob.npy"
mat3contact = add3 + \
    "general/contactmaps_numpymatrices/cacanonloc_10ang75pcntframes_.npyob.npy"

con1 = np.load(mat1contact)
con2 = np.load(mat2contact)
con3 = np.load(mat3contact)

mat1cij = add1+"general/mymatrixcij_%s.txt" % (postappend)
mat2cij = add2+"general/mymatrixcij_%s.txt" % (postappend)
mat3cij = add3+"general/mymatrixcij_%s.txt" % (postappend)

mat1 = texttomatrix(mat1cij)
mat2 = texttomatrix(mat2cij)
mat3 = texttomatrix(mat3cij)


cons_np_cij = np.zeros((len(mat1), len(mat1)))

minhas = {0: con1, 1: con2, 2: con3}
for i in range(0, len(mat1)):
    for j in range(0, len(mat1[i])):
        val1 = mat1[i][j]
        val2 = mat2[i][j]
        val3 = mat3[i][j]
        cond = abs(val1) >= cuttcij and abs(
            val2) >= cuttcij and abs(val3) >= cuttcij
        if cond:
            # avg=(val1+val2+val3)/3
            avg = (abs(val1)+abs(val2)+abs(val3))/3
            cons_np_cij[i][j] = avg

        else:
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

np.save(resfile+"_.npyob", cons_np_cij, allow_pickle=True)
with open(resfile+"_resmat", "w") as fin:
    firststring = ''
    for i in range(0, len(cons_np_cij)):
        firststring += '"V%s",' % (i+1)
    firststring = firststring[:-1]
    fin.write("%s\n" % firststring)
    for ind, val in enumerate(cons_np_cij):
        linstr = ''
        for j in val:
            linstr += "%s," % j
        fin.write('"%s",%s\n' % (ind+1, linstr[:-1]))


# below for CACAnonloc contacts, ignore above,, (should delet above)
# python jbc2016approachcij.py ../../derived/wt_r1/ ../../derived/wt_r2/ ../../derived/wt_r3/ 0.6 CA ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.6cij/cijconsensus
# python jbc2016approachcij.py ../../derived/y321a_r4/ ../../derived/y321a_r2/ ../../derived/y321a_r3/ 0.6 CA ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.6cij/cijconsensus

# python jbc2016approachcij.py ../../derived/wt_r1/ ../../derived/wt_r2/ ../../derived/wt_r3/ 0.5 CA ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/cijconsensus
# python jbc2016approachcij.py ../../derived/y321a_r4/ ../../derived/y321a_r2/ ../../derived/y321a_r3/ 0.5 CA ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.5cij/cijconsensus

# python jbc2016approachcij.py ../../derived/wt_r1/ ../../derived/wt_r2/ ../../derived/wt_r3/ 0.7 CA ../../derived/enswt/general/netob_CA/jbc2016approach_nonlocalcaca0.7cij/cijconsensus
# python jbc2016approachcij.py ../../derived/y321a_r4/ ../../derived/y321a_r2/ ../../derived/y321a_r3/ 0.7 CA ../../derived/ensy321a/general/netob_CA/jbc2016approach_nonlocalcaca0.7cij/cijconsensus


# #run netobcreator_jbc.R, and represnetations

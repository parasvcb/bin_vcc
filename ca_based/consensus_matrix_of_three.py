import numpy as np
import sys

if len(sys.argv)!=5:
    print ("Please enter correct cmd matricnumpy1 matricnumpy2 matrixnumpy3 consensusmatrixname")
    sys.exit()
prog,mat1,mat2,mat3,resfile=sys.argv

mat1=np.load(mat1)
mat2=np.load(mat2)
mat3=np.load(mat3)

cons_np=np.zeros((len(mat1),len(mat1))) 

for i in range(0,len(mat1)): 
     for j in range(0,len(mat1[i])): 
         if mat1[i][j]==mat2[i][j]==mat3[i][j] and mat1[i][j]>0: 
             cons_np[i][j]=1 
         else: 
             cons_np[i][j]=0 

np.save(resfile+"_.npyob", cons_np, allow_pickle=True)
with open(resfile+"_resmat", "w") as fin:
    firststring = ''
    for i in range(0, len(cons_np)):
        firststring += '"V%s",' % (i+1)
    firststring = firststring[:-1]
    fin.write("%s\n" % firststring)
    for ind, val in enumerate(cons_np):
        linstr = ''
        for j in val:
            linstr += "%s," % j
        fin.write('"%s",%s\n' % (ind+1, linstr[:-1]))

# python consensus_matrix_of_three.py ../../derived/wt_r1/general/contactmaps_numpymatrices/4.5ang90pcntframes_.npyob.npy  ../../derived/wt_r2/general/contactmaps_numpymatrices/4.5ang90pcntframes_.npyob.npy  ../../derived/wt_r3/general/contactmaps_numpymatrices/4.5ang90pcntframes_.npyob.npy  ../../derived/ensembledcna/dcna20_90_strict/WT_4.5ang90consensus
# python consensus_matrix_of_three.py ../../derived/wt_r1/general/contactmaps_numpymatrices/6ang90pcntframes_.npyob.npy  ../../derived/wt_r2/general/contactmaps_numpymatrices/6ang90pcntframes_.npyob.npy  ../../derived/wt_r3/general/contactmaps_numpymatrices/6ang90pcntframes_.npyob.npy  ../../derived/ensembledcna/dcna20_90_strict/WT_6ang90consensus
# python consensus_matrix_of_three.py ../../derived/y321a_r4/general/contactmaps_numpymatrices/4.5ang90pcntframes_.npyob.npy  ../../derived/y321a_r2/general/contactmaps_numpymatrices/4.5ang90pcntframes_.npyob.npy  ../../derived/y321a_r3/general/contactmaps_numpymatrices/4.5ang90pcntframes_.npyob.npy  ../../derived/ensembledcna/dcna20_90_strict/y321a_4.5ang90consensus
# python consensus_matrix_of_three.py ../../derived/y321a_r4/general/contactmaps_numpymatrices/6ang90pcntframes_.npyob.npy  ../../derived/y321a_r2/general/contactmaps_numpymatrices/6ang90pcntframes_.npyob.npy  ../../derived/y321a_r3/general/contactmaps_numpymatrices/6ang90pcntframes_.npyob.npy  ../../derived/ensembledcna/dcna20_90_strict/y321a_6ang90consensus

#go to consensus_prob.py
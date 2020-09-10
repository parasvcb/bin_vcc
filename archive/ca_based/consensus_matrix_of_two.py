import numpy as np
import sys

if len(sys.argv)!=4:
    print ("Please enter correct cmd matricnumpy1 matricnumpy2 consensusmatrixname")
    sys.exit()
prog,wt45,mt45,resfile=sys.argv

wt45=np.load(wt45)
mt45=np.load(mt45)


cons_np=np.zeros((len(wt45),len(wt45))) 

for i in range(0,len(wt45)): 
     for j in range(0,len(wt45[i])): 
         if wt45[i][j]==mt45[i][j] and wt45[i][j]>0: 
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
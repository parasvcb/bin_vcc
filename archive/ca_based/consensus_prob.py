import numpy as np,re
import sys

if len(sys.argv)!=5:
    print ("Please enter correct cmd probfile1 probfile2 probfile3 consensusmatrixname")
    sys.exit()
prog,prob1,prob2,prob3,resfile=sys.argv

def getprobhas(filename):
    with open (filename) as fin:
        dat=[i for i in fin.read().split('\n')[1:] if len(i) > 0]
    has={}
    for i in dat:
        ele=i.split('\t')
        res1,res2=list(map(int,ele[0][1:-1].split(",")))
        pair=(res1,res2) if res1<res2 else (res2,res1)
        has[pair]=float(ele[1])
    return has


prob1=getprobhas(prob1)
prob2=getprobhas(prob2)
prob3=getprobhas(prob3)

consensus_has={}
for key in prob1:
    if key in prob2 and key in prob3:
        val1=prob1[key]
        val2=prob2[key]
        val3=prob3[key]
        if 20<=val1<=90 and 20<=val2<=90 and 20<=val3<=90:
            res=sum([val1,val2,val3])/3
            consensus_has[key]=res

def writeprobfile(probhas,outfile):
    #
    keys=list(probhas.keys())
    keys.sort()
    outfile=re.sub(r'\d+pcnt.*','_prob',outfile)
    with open (outfile,"w") as fin:
        fin.write("Residpair\tP(totalframes)\n")
        for i in keys:
            if probhas[i]>0:
                fullval=probhas[i]
                fin.write("%s\t%s\n"%(i,fullval))
writeprobfile(consensus_has,resfile)

# python consensus_prob.py ../../derived/wt_r1/general/contactmaps_numpymatrices/4.5ang_prob  ../../derived/wt_r2/general/contactmaps_numpymatrices/4.5ang_prob  ../../derived/wt_r3/general/contactmaps_numpymatrices/4.5ang_prob  ../../derived/ensembledcna/dcna20_90_strict/WT_4.5ang_consensus
# python consensus_prob.py ../../derived/wt_r1/general/contactmaps_numpymatrices/6ang_prob  ../../derived/wt_r2/general/contactmaps_numpymatrices/6ang_prob  ../../derived/wt_r3/general/contactmaps_numpymatrices/6ang_prob  ../../derived/ensembledcna/dcna20_90_strict/WT_6ang_consensus
# python consensus_prob.py ../../derived/y321a_r4/general/contactmaps_numpymatrices/4.5ang_prob  ../../derived/y321a_r2/general/contactmaps_numpymatrices/4.5ang_prob  ../../derived/y321a_r3/general/contactmaps_numpymatrices/4.5ang_prob  ../../derived/ensembledcna/dcna20_90_strict/y321a_4.5ang_consensus
# python consensus_prob.py ../../derived/y321a_r4/general/contactmaps_numpymatrices/6ang_prob  ../../derived/y321a_r2/general/contactmaps_numpymatrices/6ang_prob  ../../derived/y321a_r3/general/contactmaps_numpymatrices/6ang_prob  ../../derived/ensembledcna/dcna20_90_strict/y321a_6ang_consensus
#go to comprehensive_pythonhandler_dcna.py
import sys
#from itertools import combinations
if len(sys.argv!=5):
    print ("Please Enter correct cmd 1:pdbfile 2:psffile 3:trajectory 4:resnopairs eg(135_716:...")
    sys.exit()
prog,pdb,psf,dcd,pairs=sys.argv
pairs=pairs.split(":")
pairs=[(i.split("_")[0],i.split("_")[1]) for i in pairs]
listemp=[j for i in pairs for j in i]
with open (pdb) as fin:
    datpdb=[i for i in fin.read().split("\n") if len(i)>0 and i[:4]=="ATOM"]
def getatoms(datpdb,residue):
    lis=[]
    for i in datpdb:
        atom=i[12:16]
        res=i[22:26]
        if res == residue:
            lis+=[atom]
    return lis

hasatoms={}
for i in listemp:
    hasatoms[i]=getatoms(datpdb,i)
#datafed, now compile the pairs
stringdefault='mol load psf %s dcd %s\n'%(psf,dcd)

filelist=[]
for pair in pairs:
    for atoms in hasatoms[pair[0]]:
        if atoms not in ['N','CA','C','O','HN']:
            stringdefault+='set %s [[atomselect top "resid %s and name %s"] get index]\n'%(atoms+pair[0],pair[0],atoms)
        #stringdefault+='set %s [[atomselect top "resid %s and name %s"] get index]'%(atoms,pair[0],atoms)
    for atoms in hasatoms[pair[1]]:
        if atoms not in ['N','CA','C','O','HN']:
            stringdefault+='set %s [[atomselect top "resid %s and name %s"] get index]\n'%(atoms+pair[1],pair[1],atoms)
    l1=[atoms+pair[0] for atoms in hasatoms[pair[0]]]    
    l2=[atoms+pair[1] for atoms in hasatoms[pair[1]]]
    comb = [(i,j) for i in l1 for j in l2]
    for bond in comb:
        stringdefault+='set %s [measure bond "{$%s} {$%s}" frame all]\n'%(bond[0]+'_'+bond[1],bond[0],bond[1])
        stringdefault+='set outfile [open "%s.dist" "w"]\n'%(bond[0]+'_'+bond[1])
        stringdefault+='puts $outfile %s\n'%(bond[0]+'_'+bond[1])
        filelist+=["%s.dist"%(bond[0]+'_'+bond[1])]
stringdefault+='mol delete all\nexit'

def refinetext(filename):
    with open (filename) as fin:
        dat=fin.read().split()
    with open (filename,'w') as fout:
        for ind,val in enumerate(dat):
            fout.write("%s\t%s\n"%(ind,val))


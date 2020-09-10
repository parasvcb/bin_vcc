import subprocess as sb
import os
import mdtraj as md
import sys
if len(sys.argv)!=5:
    print ("Provide correct cmd args 1. inputfile(5thcolumn as category and first as frame 1st as group) 2. WTloc 3 MTloc\
        3 outputdcdlocation")
    sys.exit()
prog,datasource,WT,Y321A,outfile=sys.argv
location={'WT':WT,'Y321A':Y321A}

def savedcd(inputfolder,lis,saveas):
    #add end in the end
    dcdfile=os.path.join(inputfolder,"productiontraj.dcd")
    pdbfile=os.path.join(inputfolder,"renumber_raw.pdb")
    dcd=md.load_dcd(dcdfile,top=pdbfile)
    t=dcd[0]
    for i in lis:
        t+=dcd[i]
    save=t[1:]
    save.save_dcd(saveas+'.dcd')
    return True

hasframes={}
with open (datasource) as fin:
    #frameno	group	variable	value	cat
    dat=[i for i in fin.read().split("\n")[1:] if len(i) >0]
    for i in dat:
        ele=i.split("\t")
        if ele[1] not in hasframes:
            hasframes[ele[1]]={}
        if ele[4] not in hasframes[ele[1]]:
            hasframes[ele[1]][ele[4]]=[]
        hasframes[ele[1]][ele[4]]+=[int(ele[0])]

for group in hasframes:
    for dcdtype in hasframes[group]:
        savedcd(location[group],hasframes[group][dcdtype],os.path.join(outfile,group+dcdtype))

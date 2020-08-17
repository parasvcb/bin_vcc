import sys
if len (sys.argv)!=5:
    print ("Please enter correct arguments: 1dcd 2pdbfile 3contacts type (0 for closest heavy atom 1 for CA) 4has_output(output)\n catdcd must be in path")
    sys.exit()
pythonprog,dcd,pdb,contype,outname=sys.argv

import pickle,subprocess,re,os,mdtraj,numpy


def contacts_computer(has,dcd,pdb):
    '''
    residue number will occasonally sgtart from 0, add 135 to them later
    has={
        (first,second):numpy.array([].dtype='float16')
    }
    '''
    print ("for dcd %s"%dcd)
    top=mdtraj.load(pdb).topology
    residues=top.n_residues
    computecon=[]
    for i in range (0,residues):
        for j in range (i,residues):
            if 0< j-i < 3:
                computecon+=[[i,j]]
    t=mdtraj.load_dcd(dcd,top=pdb)
    print ("Loaded")
    if int(contype)==0:
        print('computing closest heavy atom')
        con=mdtraj.compute_contacts(t, contacts=computecon, scheme='closest-heavy', ignore_nonprotein=True)
    else:
        print('computing CA-CA distance')
        con=mdtraj.compute_contacts(t, contacts=computecon, scheme='ca', ignore_nonprotein=True)
    print ("computed")
    object0,object1=con
    if len(has)==0:
        for i in object1:
            has[tuple(i)]=numpy.array([],dtype='float16')
    reslis=object0[0]
    for respairind in range (0,len(reslis)):
        distance_in_frames=object0[:,respairind]
        correspondingresiduepair=tuple(object1[respairind])
        has[correspondingresiduepair]=numpy.append(has[correspondingresiduepair],numpy.round(distance_in_frames*10,4))    
    return has

def trajsplitter(dcd,pdb,splitfac):
    tet=subprocess.check_output(["catdcd",dcd])
    tet=tet.decode("utf-8")
    frames=int(re.search(r'Read \d+ frames',tet).group().split()[1])
    steps=frames//splitfac
    rem=frames%splitfac
    tlis=[]
    for i in range(0,steps):
        addfac=0 if i ==0 else 1
        tlis+=[((i*splitfac)+addfac,(i+1)*splitfac)]
        i=i+1
        if rem:
            tlis+=[((i*splitfac)+addfac,(i*splitfac)+rem)]
    for ind,val in enumerate(tlis):
        createdcd=subprocess.check_output(["/home/paras/bin/./catdcd","-o","%s.tempdcd"%ind,"-first"," %s"%val[0],"-last"," %s"%val[1],dcd])

    has_mega={}
    #this should be offloaded and added if needed somewhere,a nd can be reuploaded
    for ind,val in enumerate(tlis):
        has_mega=contacts_computer(has_mega,"%s.tempdcd"%ind,pdb)
        os.remove("%s.tempdcd"%ind)
    return has_mega

has_mega=trajsplitter(dcd,pdb,1000)
append= "_ca-ca_local.hash" if int(contype) else "_closest_heavy_local.hash" 
with open (outname+append,"wb") as fin:
    pickle.dump(has_mega, fin, protocol=pickle.HIGHEST_PROTOCOL)
'''
with open (outname,"rb") as fin:
    pickle.load(fin)
'''

'''
def lframe(dcd,pdb):
    tet=subprocess.check_output(["$HOME/bin/./catdcd",dcd])
    tet=tet.decode("utf-8")
    frames=int(re.search(r'Read \d+ frames',tet).group().split()[1])
    createdcd=subprocess.check_output(["/home/paras/bin/./catdcd","-o","testpull.dcd","-first"," %s"%(frames-2),"-last"," %s"%(frames),dcd])
    t = mdtraj.load_dcd("testpull.dcd",top=pdb)
    lastframedcd=mdtraj.load_dcd("testpull.dcd",top=pdb,frame=t.n_frames-1)
    lastframedcd.remove_solvent().save_pdb("/tmp/last.pdb")
    return True
'''

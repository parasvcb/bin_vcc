import subprocess as sb
import os
import sys
if len(sys.argv)!=5:
    print ("Provide correct cmd args 1. inputfile(5thcolumn as category and first as frame 1st as group) 2. group folder WT 3 groupfolder MT\
        3 outputdcdlocation")
    sys.exit()
prog,datasource,WT,Y321A,outfile=sys.argv
location={'WT':WT,'Y321A':Y321A}

def savedcd(inputfolder,lis,saveas):
    #add end in the end
    stringinput=''
    count=1
    for i in lis:    
        with open (os.path.join('%s'%inputfolder,'%s.pdb'%i)) as fin:
            stringinput+='MODEL        %s'%count
            stringinput+='\n'+'\n'.join([i for i in fin.read().split("\n") if len(i)>0 and i[0:4]=='ATOM'])
            stringinput+='\n'+'TER    8940      ASN A 716\nENDMDL'+'\n'
            count+=1
    stringinput+='END'
    with open(saveas+".pdb",'w') as fin:
        fin.write("%s"%stringinput)
    res=sb.call(['mdconvert',saveas+".pdb",'-o',saveas+'.dcd'])
    #print (stringinput)
    #print (saveas)
    os.remove(saveas+'.pdb')
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

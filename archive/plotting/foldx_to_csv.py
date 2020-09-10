import sys
import os
if len(sys.argv)!=4:
    print ("please neter corect directories wtdir y321adir outfile")
    sys.exit()

prog,simdirwt,simdirmt,outfile=sys.argv

def has_updater(fileadd,has):
    with open (fileadd) as fin:
        dat=fin.read().split('\n')[9:-1]
    for interactionpair in dat:
        elements=interactionpair.split("\t")
        pair=[elements[1],elements[2]]
        intenergy=[elements[-3]]
        pair.sort()
        key="".join(pair)
        if key not in has:
            has[key]=[]
        has[key]+=intenergy
    return has

def filewriter(haswt,hasmt,filename):
    with open (filename,"w") as fin:
        keys=list(haswt.keys())
        keys.sort()
        fin.write("domgroup,int_energy,type\n")
        for i in keys:
            for j in haswt[i]:
                fin.write("%s,%s,WT\n"%(i,j))
            for j in hasmt[i]:
                fin.write("%s,%s,Y321A\n"%(i,j))

dom4_has_wt={}
for files in os.listdir(simdirwt):
    if files[:7] == 'Summary':
        dom4_has_wt=has_updater(os.path.join(simdirwt,files),dom4_has_wt)

dom4_has_mt={}
for files in os.listdir(simdirmt):
    if files[:7] == 'Summary':
        dom4_has_mt=has_updater(os.path.join(simdirmt,files),dom4_has_mt)

filewriter(dom4_has_wt,dom4_has_mt,outfile)

import sys,os
if len(sys.argv)!=4:
    print ("please neter corect directories wtdir y321adir outfileappend")
    sys.exit()

'''
ATOM     20  CA  THR C 136      27.151  83.185   8.041  1.00  0.00
ATOM     34  CA  ASN C 137      27.058  85.741   5.274  2.00  0.00
ATOM     48  CA  THR C 138      30.856  86.580   5.238  3.00  0.00
'''
prog,wtdir,mtdir,outappend=sys.argv
def fileparser(inpdir):
    print (inpdir)
    has={}
    fcount=[]
    for fil in [i for i in os.listdir(inpdir) if i.split(".")[-1]=='pdb']:
        #print (fil)
        fcount+=[fil.split("_")[0]]
        #fileframenoumber
        with open (inpdir+fil) as fin:
            dat=[i for i in fin.read().split('\n') if len(i)>0 and i[:4]=="ATOM"]
            lis=[]
        for line in dat:
            ele=line.split()
            resid=ele[5]
            chain=ele[4]
            resname=ele[3]
            lis+=[(chain,chain+'_'+resid+'_'+resname)]
            #chain, chainID+residueID+resname
        chaing=[i[0] for i in lis]
        chaing=list(set(chaing))
        chaing.sort()
        key=chaing[0]+chaing[1]
        if key not in has:
            has[key]=[[],{}]
        has[key][0]+=[len(lis)]
        for ele in lis:
            if ele[1] not in has[key][1]:
                has[key][1][ele[1]]=0
            has[key][1][ele[1]]+=1
    '''
    has[domgroup]=[[listofnoofiresiduesinfrmaes],{chain_residue_resid:count(depictingnumber_ofcocurences)}]
    '''
    return has,len(set(fcount))

wthas,fcount=fileparser(wtdir)
mthas,fcount=fileparser(mtdir)
print(fcount)
#print (wthas)


#sed -e  's/^/source\/ensy321a\/dcd_to_pdb\//' -i ensy321apdblis

def writetofile(haswt,hasmt,append):
    keys=list(haswt.keys())
    keys.sort()
    with open(append+"totalpairs.csv",'w') as fin:
        fin.write("domgroup,rescount,type\n")
        for domgroup in keys:
            for noofintresinframes in haswt[domgroup][0]:
                fin.write("%s,%s,WT\n"%(domgroup,noofintresinframes))
            for noofintresinframes in hasmt[domgroup][0]:
                fin.write("%s,%s,MT\n"%(domgroup,noofintresinframes))
    
    with open(append+"reswise.csv",'w') as fin:
        fin.write("domgroup,resdet,frequency,type\n")
        for domgroup in keys:
            resids=list(set(list(haswt[domgroup][1].keys())+list(hasmt[domgroup][1].keys())))
            resids.sort()
            #print (resids)
            for res in resids:
                wtfreq=haswt[domgroup][1][res]/fcount if res in haswt[domgroup][1] else 0
                mtfreq=hasmt[domgroup][1][res]/fcount if res in hasmt[domgroup][1] else 0
                if abs(wtfreq-mtfreq)>0.25:
                    print (wtfreq,mtfreq,wtfreq-mtfreq)    
                    fin.write("%s,%s,%s,WT\n"%(domgroup,res,wtfreq))
                    fin.write("%s,%s,%s,MT\n"%(domgroup,res,mtfreq))
    
writetofile(wthas,mthas,outappend)


'''
df=read.csv('ialignresreswise.csv')
gg <- ggplot(data=df, aes(x=resdet,y=frequency,fill=type), size=0.5,alpha=0.7)
gg <- gg+geom_bar(stat='identity',position=position_dodge())
gg <- gg + facet_wrap(~domgroup,ncol=1,scales='free_x')
gg <- gg + labs(y="Freq", x="Residue")
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 8, angle = 45), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))

'''
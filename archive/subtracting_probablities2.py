import sys
import os
#from matplotlib_venn import venn3,venn2
#import matplotlib.pyplot as plt
if len(sys.argv)!=3:
    print ("Please enter correct cmd arguments probfile, outputfile")
    sys.exit()
prog,probfile,outfile=sys.argv

def has_updater(fileadd,has):
    with open (fileadd) as fin:
        dat=[i for i in fin.read().split('\n')[1:] if len(i)>0]
    for interactionpair in dat:
        elements=interactionpair.split("\t")
        diff=float(elements[-1])
        r1,r2=map(int,elements[0][1:-1].split(",")) 
        intprobability=float(elements[1])
        pair=[r1+135,r2+135]
        pair.sort()
        has[tuple(pair)]=[diff,intprobability]
        #adding difference interactions in first and last halves as 0th elemnet
    return has

def filewriter_csv(has,filename1):
    superkeys=[]
    with open (filename1,"w") as fin:
        #fin.write("Group,index,difference\n")
        fin.write("less30\tless25\tless20\tless15\tless10\n")
        values=[i[0] for i in list(has.values())]
        less20=[i for i in values if abs(i)<=20]
        #that means value difference should be in between -20 to 20
        less15=[i for i in values if abs(i)<=15]
        less10=[i for i in values if abs(i)<=10]
        less25=[i for i in values if abs(i)<=25]
        less30=[i for i in values if abs(i)<=30]
        fin.write("%s\t%s\t%s\t%s\t%s\n"%(len(less30)/len(values),len(less25)/len(values),len(less20)/len(values),len(less15)/len(values),len(less10)/len(values)))
        fin.write("\n")
    #this function will write the contents and pairs in threshold for all the replicates

has_interactions={}
has_interactions=has_updater(probfile,has_interactions)
filewriter_csv(has_interactions,outfile+"containment_in_thresholdbarriers.csv")


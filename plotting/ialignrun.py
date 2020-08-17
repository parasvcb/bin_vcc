import subprocess
import sys,os
if len(sys.argv)!=5:
    print ("Please enter correct cmd 1.ialign execyutable 2. inputdir 3outputdiri 4outfile")
    sys.exit()
prog,src,inp,out,of=sys.argv
counter=0
with open(of,'w') as fin:
    for i in os.listdir(inp):
        counter+=1
        inpt='ensy321a/'
        out='ensy321aout/'
        if counter<=119:
            fin.write("%s -w %s %s %s -a 0 -g &\n"%(src,out,inpt+i,inpt+i))
        else:
            fin.write("%s -w %s %s %s -a 0 -g\n"%(src,out,inpt+i,inpt+i))
        if counter==120:
            counter=1
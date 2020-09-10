import sys
import os
from progress.bar import Bar
if len(sys.argv) != 3:
    print("Please enter correct args inpfile, outfile willbenamedas 5domain.pdb")
    sys.exit()
prog, inputfile, outfile_prefix = sys.argv
has_5dom = {}


for i in range(21, 717):
    if 135 <= i < 277:
        has_5dom[i] = 'C'
    elif 325 <= i < 455:
        has_5dom[i] = 'C'
    elif 455 <= i < 582:
        has_5dom[i] = 'H'
    elif 277 <= i < 325:
        has_5dom[i] = 'M'
    elif 582 <= i < 717:
        has_5dom[i] = 'P'
    else:
        has_5dom[i] = 'U'

def pdb_renumber(source, outfile):
    # pdb 22 to 25,
    newdat5 = []
    for i in source:
        if len(i) > 10 and i[:4] == 'ATOM':
            newdat5 += [i[:21]+"%s" % (has_5dom[int(i[22:26])])+i[22:]]
        else:
            newdat5 += [i]
    with open(outfile+'_4dom.pdb', "w") as fin:
        for i in newdat5:
            fin.write("%s\n" % i)

with open(inputfile) as fin:
    dat = fin.read().split("\n")
pdb_renumber(dat, outfile_prefix)

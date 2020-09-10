import sys
import os
from progress.bar import Bar
if len(sys.argv) != 3:
    print("Please enter correct args pdbdirectory, outputdirectory (willbenamedas 5domain.pdb and 3domain.pdb ahead of input name)")
    sys.exit()
prog, inputfile, outfile_prefix = sys.argv
has_5dom = {}
has_3dom = {}

for i in range(21, 717):
    if 135 <= i < 277:
        has_5dom[i] = 'C'
        has_3dom[i] = 'C'
    elif 325 <= i < 455:
        has_5dom[i] = 'C'
        has_3dom[i] = 'C'
    elif 455 <= i < 582:
        has_5dom[i] = 'H'
        has_3dom[i] = 'H'
    elif 277 <= i < 325:
        has_5dom[i] = 'M'
        has_3dom[i] = 'C'
    elif 582 <= i < 717:
        has_5dom[i] = 'P'
        has_3dom[i] = 'P'
    else:
        has_5dom[i] = 'U'
        has_3dom[i] = 'U'


def pdb_renumber(source, outfile):
    # pdb 22 to 25,
    newdat5 = []
    newdat3 = []
    for i in source:
        if len(i) > 10 and i[:4] == 'ATOM':
            newdat5 += [i[:21]+"%s" % (has_5dom[int(i[22:26])])+i[22:]]
            newdat3 += [i[:21]+"%s" % (has_3dom[int(i[22:26])])+i[22:]]
        else:
            newdat5 += [i]
            newdat3 += [i]
    with open(outfile+'_4dom.pdb', "w") as fin:
        for i in newdat5:
            fin.write("%s\n" % i)

inputlist = [i for i in os.listdir(os.path.join(inputfile)) if i.split(".")[-1]=='pdb']
bar = Bar('Processing', max=len(inputlist))
for i in inputlist:
    with open(os.path.join(inputfile, i)) as fin:
        dat = fin.read().split("\n")
    pdb_renumber(dat, os.path.join(outfile_prefix)+i.split('.')[0])
    bar.next()
bar.finish()

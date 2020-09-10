import re
import sys
if len(sys.argv) != 4:
    print("Please enter correct cmd args 1.cijwt 2cijmt 3.residuecomfile 3outputfile")
    sys.exit()
prog, filenamesource, rescom, out = sys.argv

cijval = {}
with open(filenamesource, 'r') as fin:
    for line in [i for i in fin.read().split("\n")[1:] if len(i) > 0]:
        ele = line.split("\t")
        pair = list(map(int, ele[:2]))
        cjwt, cjmt, cjsub = list(map(float, ele[2:]))
        cijval[tuple(pair)] = [cjwt, cjmt, cjsub]


def vmd_displayer(centerfile, out):
    with open(centerfile) as fin:
        stringfile = re.sub(r'{', '', fin.read())
        stringfile = re.sub(r'}', '', stringfile)
        dat = [i for i in stringfile.split("\n") if len(i) > 0]
    hascom = {int(i.split()[0]): " ".join(i.split()[1:]) for i in dat}
    with open(out, 'w') as vmd:
        for pair in cijval:
            wtcoupling, mtcoupling, subcoupling = cijval[pair]
            if (abs(wtcoupling) >= 0.4 or abs(mtcoupling) >= 0.4) and abs(subcoupling) >= 0.4:
                # means correlations should be stronger in either wt[0] and y321a[1] and their difference[2] should also be greater than 0.4
                # red ccolor for gain in correlation
                # pink for gain in anticorrelation
                # blue for loss in correlation
                # green for loss in anticorrelation
                # color = 'red' if mtcoupling>=0.4 else \
                #     'pink' if wtcoupling>-0.4 and mtcoupling<=-0.4 else \
                #     'blue' if wtcoupling>=0.4 and mtcoupling<=-0.4 else
                res1, res2 = pair
                value = cijval[pair][2]
                color = 'blue' if value < 0 else 'red'
                cylinderwidth = abs(value)*1
                # above value can be varied to change width for substraction couplings
                string = '''
                draw color %s
                draw sphere {%s} radius 0.5 resolution 25
                draw color %s
                draw sphere {%s} radius 0.5 resolution 25
                draw color %s
                draw cylinder {%s} {%s} radius %s resolution 25
                ''' % (color, hascom[res1], color, hascom[res2], color, hascom[res1], hascom[res2], cylinderwidth)
                vmd.write("%s\n" % string)


vmd_displayer(rescom, out)

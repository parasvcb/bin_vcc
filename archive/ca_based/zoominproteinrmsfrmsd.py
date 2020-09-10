import sys
if len(sys.argv) != 2:
    print("Please enter correct cmd args 1 inputproteinfile")
    sys.exit()
prog, inp = sys.argv
'''
_all.pdf                _dom2_1_rmsd.txt        _dom2_2_rmsf.txt        _dom2.pdf               _dom3_rmsd.txt          _dom4_rmsf.txt          rmsd_all.pdf
_all_rmsd.txt           _dom2_1_rmsf.txt        _dom2_3.pdf             _dom2_rmsd.txt          _dom3_rmsf.txt          _membindmotif.pdf       rmsf_6domains.pdf
_all_rmsf.txt           _dom2_2.pdf             _dom2_3_rmsd.txt        _dom2_rmsf.txt          _dom4.pdf               _membindmotif_rmsd.txt  
_dom2_1.pdf             _dom2_2_rmsd.txt        _dom2_3_rmsf.txt        _dom3.pdf               _dom4_rmsd.txt          _membindmotif_rmsf.txt  
'''
dom2_1lis = list(range(135, 277))
dom2_2lis = list(range(325, 455))
dom2_3lis = list(range(455, 582))
dom2lis = dom2_1lis+dom2_2lis+dom2_3lis
dom3lis = list(range(277, 325))
dom4lis = list(range(582, 717))
memlis = list(range(292, 311))
alllis = list(range(135, 717))
dom2_1list = []
dom2_2list = []
dom2_3list = []
dom2list = []
dom3list = []
dom4list = []
memlist = []
alllist = []


def liswriter(lis, filename):
    filename = "/".join(inp.split("/")[:-1])+"/_"+filename+"_rmsfzoom.txt"
    with open(filename, 'w') as fin:
        fin.write("rmsf\tresid\tgrp\tprotgrp\n")
        for i in lis:
            fin.write("%s\t%s\t%s\t%s\n" % (i[0], i[1], i[2], i[3]))
    return


with open(inp) as fin:
    dat = [i for i in fin.read().split("\n")[1:] if len(i) > 0]
    for i in dat:
        rmsf, resid, grp, protgrp = i.split("\t")
        resid = int(resid)
        if resid in dom2_1lis:
            dom2_1list += [[rmsf, resid, grp, "dom2_1"]]
        if resid in dom2_2lis:
            dom2_2list += [[rmsf, resid, grp, "dom2_2"]]
        if resid in dom2_3lis:
            dom2_3list += [[rmsf, resid, grp, "dom2_3"]]
        if resid in dom2lis:
            dom2list += [[rmsf, resid, grp, "dom2"]]
        if resid in dom3lis:
            dom3list += [[rmsf, resid, grp, "dom3"]]
        if resid in dom4lis:
            dom4list += [[rmsf, resid, grp, "dom4"]]
        if resid in memlis:
            memlist += [[rmsf, resid, grp, "membindmotif"]]
        if resid in alllis:
            alllist += [[rmsf, resid, grp, "whole_protein"]]
liswriter(dom2_1list, "dom2_1")
liswriter(dom2_2list, "dom2_2")
liswriter(dom2_3list, "dom2_3")
liswriter(dom3list, "dom3")
liswriter(dom2list, "dom2")
liswriter(dom4list, "dom4")
liswriter(memlist, "membindmotif")
liswriter(alllist, "all")


# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/wt_r1/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/wt_r2/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/wt_r3/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/y321a_r1/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/y321a_r2/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/y321a_r3/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/enswt23repgen/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/enswtgen/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/ensy321a23repgen/_all_rmsf.txt
# python bin/zoominproteinrmsfrmsd.py derived/no_pro_simulation/ensy321agen/_all_rmsf.txt

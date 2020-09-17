# This folder will tease out the general layouts and details of the network methodology and will specfilaly highlight
# 1. Distributions and their density maps/histograms of cij values, (6 in one panel).
# 2. Number of pairs with cij >=0.6 in all the replicates (6 values in return).
# 3. Among them how many were from consensus and how many were distance based (i think i know them) and their fraction from total pairs)
# 4. Get a common set of wildtype and mutant and subtract their values from wildtype and display distributions.
# 5. focus on fraction that will make a difference in general (non zero likely)

import sys
import os
import re
import numpy as np
if len(sys.argv) != 4:
    print("please enter correct cmd args 1wt dir 2 y321adir 3outdir")
    sys.exit()
prog, wt, yt, outdir = sys.argv
wtdir = ['WT_r1', 'WT_r2', 'WT_r3']
ytdir = ['Y321A_r4as1', 'Y321A_r2', 'Y321A_r3']
exDir = 'output_process_renumberpdb_aligned'


def texttomatrix(matfile):
    with open(matfile) as fin:
        dat = fin.read()
    dat = re.sub(r'"', '', dat)
    dat = re.sub(r'V', '', dat)
    with open("testmat", 'w') as fin:
        fin.write("%s" % dat)
    a = np.loadtxt(open("testmat"), delimiter="\t", skiprows=1)
    a = np.array([i[1:] for i in a])
    a = a[np.triu_indices(len(a), k=1)]
    # get only upper half of matrix without 1 diagnoal values
    os.remove('testmat')
    return np.abs(a)


def cijmatrixtodf_and_obs2(has, outfile):
    with open(outfile, 'w') as fout:
        fout.write("group,cij\n")
        for i in has:
            for j in has[i]:
                fout.write("%s,%s\n" % (i, round(j, 3)))
    for i in has:
        countgrt1 = len(np.where(has[i] >= 0.1)[0])
        countgrt6 = len(np.where(has[i] >= 0.6)[0])
        print("for group %s, >=0.6 values were:%s (%s freq from values gretar than 0.1 (%s pairs) and %s (includin all pairs, which were %s))"
              % (i, countgrt6, round(countgrt6/countgrt1, 3), countgrt1, round(countgrt6/len(has[i]), 3), len(has[i])))


def mattohas(mat):
    def valgrtreturn(val, denom):
        if val:
            return round(val/denom, 3)
        else:
            return 0
    has = {}
    for i in range(0, len(mat)):
        for j in range(0, len(mat[i])):
            if i != j:  # removing diagnoal elements
                key = [i, j]
                key.sort()  # removing repetitions
                has[tuple(key)] = abs(mat[i, j])
    values = np.asarray(list(has.values()))
    totalpairs = len(values)
    print(totalpairs)
    valgrt1 = len(np.where(values >= 0.1)[0])
    valgrt2 = len(np.where(values >= 0.2)[0])
    valgrt3 = len(np.where(values >= 0.3)[0])
    valgrt4 = len(np.where(values >= 0.4)[0])
    valgrt5 = len(np.where(values >= 0.5)[0])
    valgrt6 = len(np.where(values >= 0.6)[0])
    valgrt7 = len(np.where(values >= 0.7)[0])
    valgrt8 = len(np.where(values >= 0.8)[0])
    valgrt9 = len(np.where(values >= 0.9)[0])
    print('in mattohas()')
    string = '''Total greater than 0.1 were: %s (%s)
    Total greater than 0.2 were: %s (%s)
    Total greater than 0.3 were: %s (%s)
    Total greater than 0.4 were: %s (%s)
    Total greater than 0.5 were: %s (%s)
    Total greater than 0.6 were: %s (%s)
    Total greater than 0.7 were: %s (%s)
    Total greater than 0.8 were: %s (%s)
    Total greater than 0.9 were: %s (%s)
    ''' % (valgrt1, valgrtreturn(valgrt1, totalpairs), valgrt2, valgrtreturn(valgrt2, totalpairs),
           valgrt3, valgrtreturn(valgrt3, totalpairs), valgrt4, valgrtreturn(
               valgrt4, totalpairs),
           valgrt5, valgrtreturn(valgrt5, totalpairs), valgrt6, valgrtreturn(
               valgrt6, totalpairs),
           valgrt7, valgrtreturn(valgrt7, totalpairs), valgrt8, valgrtreturn(
               valgrt8, totalpairs),
           valgrt9, valgrtreturn(valgrt9, totalpairs))
    print(string)
    return {i: has[i] for i in has if has[i] >= 0.6}


w1 = texttomatrix(os.path.join(
    wt, exDir, wtdir[0], 'sim_data', 'mymatrixcij_CA.txt'))
w2 = texttomatrix(os.path.join(
    wt, exDir, wtdir[1], 'sim_data', 'mymatrixcij_CA.txt'))
w3 = texttomatrix(os.path.join(
    wt, exDir, wtdir[2], 'sim_data', 'mymatrixcij_CA.txt'))
y1 = texttomatrix(os.path.join(
    yt, exDir, ytdir[0], 'sim_data', 'mymatrixcij_CA.txt'))
y2 = texttomatrix(os.path.join(
    yt, exDir, ytdir[1], 'sim_data', 'mymatrixcij_CA.txt'))
y3 = texttomatrix(os.path.join(
    yt, exDir, ytdir[2], 'sim_data', 'mymatrixcij_CA.txt'))
has = {'WT_r1': w1, 'WT_r2': w2, 'WT_r3': w3,
       'Y321A_r1': y1, 'Y321A_r2': y2, 'Y321A_r3': y3}
# mymatrixcij_CA.txt
# CONSENSUS_Cij_MATRIX_.npyob.npy
cijmatrixtodf_and_obs2(has, os.path.join(outdir, 'raw_distribution.csv'))

print('WTcons')
wcons = mattohas(np.load(os.path.join(
    wt, exDir, 'CONSENSUS_Cij_MATRIX_.npyob.npy')))

print('Y321Acons')
ycons = mattohas(np.load(os.path.join(
    yt, exDir, 'CONSENSUS_Cij_MATRIX_.npyob.npy')))

commonkeys = set(list(wcons.keys())+list(ycons.keys()))
print(list(ycons.keys())[:5])
diffval = []
with open(outdir+'difference_finaledges.csv', 'w') as fout:
    fout.write('values\n')
    for i in commonkeys:
        wval = wcons[i] if i in wcons else 0
        yval = ycons[i] if i in ycons else 0
        diff = wval-yval
        diffval += [diff]
        fout.write('%s\n' % round(diff, 3))
print("For final consensus")
print("For Y321A -> Pairs were %s:" % (len(ycons.keys())))
print("For WT -> Pairs were %s:" % (len(wcons.keys())))
print("common set -> Pairs were %s:" % (len(commonkeys)))
vals, bins = np.histogram(
    diffval, bins=list(np.arange(-1.1, 1.1, 0.1)))

for ind, val in enumerate(bins):
    if ind < len(bins)-1:
        print("For the bin ranging %s to %s, values were %s(%s)"
              % (round(bins[ind], 2), round(bins[ind+1], 2), vals[ind], np.round(vals[ind]/sum(vals), 2)))
# python calculations.py ../../../../WT_series/ ../../../../Y321A_series/ ../../../../plotting/ >calculationlog.log
# divide has_out with nonzero values of has, to get 2 in freq

# 3rd i do have

# for 4th read the consensus cij matrix and compare the non zero values.
# add the to collective hash and plot them,
# >>> 4960/((109049+111276+112715)/3)*100
# 4.467931779966371
# >>> 5294/((113891+120488+110337)/3)
# 0.04607270912867404

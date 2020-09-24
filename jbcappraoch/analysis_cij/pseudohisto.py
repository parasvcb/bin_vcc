with open('../../../../plotting/raw_distribution.csv') as fin:
    dat = [i for i in fin.read().split('\n')[1:] if len(i) > 0]

binshas = {(0, 0.1): [], (0.1, 0.2): [], (0.2, 0.3): [], (0.3, 0.4): [], (0.4, 0.5): [
], (0.5, 0.6): [], (0.6, 0.7): [], (0.7, 0.8): [], (0.9, 1.0): [], (1.0, 1.1): []}
hascustom = {'Y321A_r1': {}, 'Y321A_r2': {},
             'Y321A_r3': {}, 'WT_r1': {}, 'WT_r2': {}, 'WT_r3': {}}
hascustom['Y321A_r1'] = {i: 0 for i in binshas}
hascustom['Y321A_r2'] = {i: 0 for i in binshas}
hascustom['Y321A_r3'] = {i: 0 for i in binshas}
hascustom['WT_r1'] = {i: 0 for i in binshas}
hascustom['WT_r2'] = {i: 0 for i in binshas}
hascustom['WT_r3'] = {i: 0 for i in binshas}
# print(dat[:5])
for record in dat:
    ele = record.split(',')
    cij = float(ele[1])
    key = ele[0]
    for j in binshas:
        if j[0] <= cij < j[1]:
            hascustom[key][j] += 1
            break
keys = list(binshas.keys())
keys.sort()
with open('../../../../plotting/raw_distribution_pythonpseudoHist_with0bin.csv', 'w') as fout:
    fout.write('bins,reptype,values,frequency\n')
    for reptype in hascustom:
        repsum = sum(hascustom[reptype].values())
        print(reptype, repsum)
        for bins in keys:
            binval = hascustom[reptype][bins]
            fout.write("%s,%s,%s,%s\n" %
                       (bins[0], reptype, binval, round(binval / repsum, 2)))

keys = keys[1:]
print('newhascutsomwithout0bin')
newhascustom = {}
for reptype in hascustom:
    if reptype not in newhascustom:
        newhascustom[reptype] = {}
    for bins in hascustom[reptype]:
        if bins != (0, 0.1):
            newhascustom[reptype][bins] = hascustom[reptype][bins]
with open('../../../../plotting/raw_distribution_pythonpseudoHist_without0bin.csv', 'w') as fout:
    fout.write('bins,reptype,values,frequency\n')
    for reptype in newhascustom:
        repsum = sum(newhascustom[reptype].values())
        print(reptype, repsum)
        for bins in keys:
            binval = newhascustom[reptype][bins]
            fout.write("%s,%s,%s,%s\n" %
                       (bins[0], reptype, binval, round(binval / repsum, 2)))

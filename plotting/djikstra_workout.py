import networkx as nx
import numpy as np
import sys
if len(sys.argv) != 4:
    print("please 4 args 1.adj matrix WT 2. adjmatrix Y321A 2. outputappend")
    sys.exit()
prog, matfileWT, matfileY321A, outappend = sys.argv
'''
The data tsructure i will be needing for this:
1. matrix excel, two hashes and analysis,
	both the groups, if node pairs, (source dest), then path and weight else false and 0
2. For R plots,
	path length distributions
		group1 and group2 in one column and WT and Y321A in one column, followed by shortes path length values
	ndoe degenrecay values
		group1 and group2 in one column and WT and Y321A in one column, nodes and theri value sin another clumns
		sort according to WT,
		sort by bar chart
'''
wt = np.loadtxt(matfileWT)
y321a = np.loadtxt(matfileY321A)
wtg = nx.from_numpy_matrix(wt, create_using=None)
mtg = nx.from_numpy_matrix(y321a, create_using=None)
sourcegroup = {'Tyr321': [321]}
# 'loop17': list(range(361, 370)),  # yellow
# 'loop8': list(range(233, 240)),  # red
# 'loop21': [421, 422, 423, 424],  # orange
# 'allSource': list(range(361, 370)) + list(range(233, 240)) + [421, 422, 423, 424]}

destinationgroup = {'loop38': list(range(578, 586)),
                    'loop11': list(range(292, 311)),
                    'prismdomain': list(range(586, 717)),
                    'val257': [257],
                    'allSinkwoprism': list(range(578, 586)) + list(range(292, 311))+[257],
                    'allSink': list(range(578, 586)) + list(range(292, 311))+list(range(586, 717))+[257]}


def compute_distances(g, lis1, lis2):
        # NS, SOURCE LIS, DEST LIST
    residfac = 135
    has = {}
    for i in lis1:
        for j in lis2:
            source = i-residfac
            dest = j-residfac
            try:
                length = nx.dijkstra_path_length(
                    g, source=source, target=dest, weight='weight')
                path = nx.dijkstra_path(
                    g, source=source, target=135, weight='weight')
                has[(source+residfac, dest+residfac)
                    ] = [length, [k+135 for k in path]]
            except Exception as E:
                has[(source+residfac, dest+residfac)] = [0, ['False']]
    return has


def nodedegeneracy(list, key):
    has = {}
    src, dest = key.split("_")
    exlcude = [sourcegroup[src]+destinationgroup[dest]]
    # exlcude source and destinations groups
    all = [j for i in list for j in i if j not in exlcude]
    has = {i: [all.count(i), len(list), all.count(i)/len(list)]
           for i in set(all)}
    # has will have node as key, its count, total paths, and frequency
    return has


has_wt = {}
has_mt = {}
for i in sourcegroup:
    for j in destinationgroup:
        new_key = i+"_"+j  # source and destinations
        has_wt[new_key] = compute_distances(
            wtg, sourcegroup[i], destinationgroup[j])
        has_mt[new_key] = compute_distances(
            mtg, sourcegroup[i], destinationgroup[j])

# hashes has been filled
# Give three tsv's
comp = open(outappend+"comprehensive_paths.tsv", 'w')  # indivisual group name
comp.write("g_Src\tg_Dest\tsource\tsink\tpath\tlength\ttype\n")

pathdist = open(outappend+"plot_paths_distrib.tsv", 'w')  # common group name
pathdist.write("g_Src_Dest\tsource\tsink\tlength\tpaths\ttype\n")

nodedist = open(outappend+"plot_node_distrib.tsv", 'w')  # common group name
nodedist.write("g_Src_Dest\tpaths\tnode\tnodeval\tnodefreq\ttag\ttype\n")

for group in has_wt:
    sg, dg = group.split("_")
    wtpaths = []
    mtpaths = []
    pathcwt = [1 for ti in has_wt[group] if has_wt[group]
               [ti][1][0] != 'False']  # ti is resid pair, tj is val
    pathcmt = [1 for ti in has_mt[group] if has_mt[group]
               [ti][1][0] != 'False']  # ti is resid pair, tj is val
    # print(sum(pathcwt), len(has_wt[group]), sum(pathcmt), len(has_mt[group]))
    for pairs in has_wt[group]:
        lengthwt, pathwt = has_wt[group][pairs]
        lengthmt, pathmt = has_mt[group][pairs]
        comp.write("%s\t%s\t%s\t%s\t%s\t%s\tWT\n" %
                   (sg, dg, pairs[0], pairs[1], pathwt, lengthwt))
        comp.write("%s\t%s\t%s\t%s\t%s\t%s\tY321A\n" %
                   (sg, dg, pairs[0], pairs[1], pathmt, lengthmt))

        if pathwt[0] != "False":
            pathdist.write("%s\t%s\t%s\t%s\t%s\tWT\n" %
                           (group, pairs[0], pairs[1], lengthwt, sum(pathcwt)))
        if pathmt[0] != "False":
            pathdist.write("%s\t%s\t%s\t%s\t%s\tY321A\n" %
                           (group, pairs[0], pairs[1], lengthmt, sum(pathcmt)))
        # add up all the path and send for verification later
            wtpaths += [pathwt]
            mtpaths += [pathmt]
    wtfreq = nodedegeneracy(wtpaths, group)
    mtfreq = nodedegeneracy(mtpaths, group)
    common_nodes = set(list(wtfreq.keys())+list(mtfreq.keys()))
    for nodes in common_nodes:
        wtnodec, wtpathc, wtnodef = wtfreq[nodes] if nodes in wtfreq else [
            0, 0, 0]
        mtnodec, mtpathc, mtnodef = mtfreq[nodes] if nodes in mtfreq else [
            0, 0, 0]
        tag = 'plot' if (wtnodef > 0.1 and mtnodef > 0.1) or abs(
            wtnodef-mtnodef) >= 0.25 else 'data'
        nodedist.write("%s\t%s\tres_%s\t%s\t%s\t%s\tWT\n" %
                       (group, wtpathc, nodes, wtnodec, wtnodef, tag))
        nodedist.write("%s\t%s\tres_%s\t%s\t%s\t%s\tY321A\n" %
                       (group, mtpathc, nodes, mtnodec, mtnodef, tag))


# python djikstra_workout.py ../../derived/cnapathobjects/adj_jbc_wt.txt ../../derived/cnapathobjects/adj_jbc_y321a.txt jbc_shortestpaths_ &
# python djikstra_workout.py ../../derived/cnapathobjects/adj_cij0.6_4.5_50_wt.txt ../../derived/cnapathobjects/adj_cij0.6_4.5_50_y321a.txt 4.5_50_cij0.6_shortestpaths_ &

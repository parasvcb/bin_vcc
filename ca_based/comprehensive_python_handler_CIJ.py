import os
from shutil import rmtree
import sys
if len(sys.argv) != 9:
    print(len(sys.argv))
    print("""
    Please enter the correct command line arguements:
    1: ensemble/replicates, if replicates then also mention frame count with 10000:20000 or more (eg replicates_20000:30000)
    2:dcd 3:pdb 4:outputdir(with sub classification of netob general analysis)
    5: ultrsmalldcd 6:numpycontacts 7:aftreshstart(0 for resume with preexisting files, 1 delete all and proceed,
    8: cores)
    """)
    sys.exit()


def createdir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


prog,systemtype, dcd, pdb, mainout, ultrsmalldcd, numpyhas, afreshstart, cores = sys.argv
replicates=True if "rep" in systemtype else False
if replicates:
    frames=systemtype.split('_')[1]
if int(afreshstart):
    if os.path.isdir(mainout):
        rmtree(mainout)
        createdir(mainout)
dcddir="/".join(dcd.split("/")[:-1])

genout=mainout+"general/"
numpymat=genout+"contactmaps_numpymatrices/"
netobout=genout+"netob_CA/"
netobanaout=genout+"netobanalysis_CA/"

createdir(genout)
createdir(numpymat)
createdir(netobout)

import subprocess
import multiprocessing
cijob=genout+"CIJ_CA.ob"
xyzob=genout+"XYZob_CA.ob"

def numpycmaprunner(ultrsmalldcd, cijob, numpymat, outname):
    subprocess.check_output(
        ["Rscript", "cmap_dummyob_creator_s2_1.R", ultrsmalldcd, cijob, outname, numpymat])

def netobcreator(cijob, cmapob, cijcutt, outname):
    subprocess.check_output(
        ["Rscript", "netob_creator_with_cmap_s3.R", cijob, cmapob, cijcutt, outname])

def netobanalysis(netob, outname, pdb, cijob):
    subprocess.check_output(
        ["Rscript", "netob_analysis_s4.R", netob, outname, pdb, cijob])

if not (os.path.isfile(cijob) and os.path.isfile(xyzob)):
    step1=subprocess.check_output(
        ["Rscript", "trajectory_general_analysis_s1.R", dcd, pdb, genout, cores])
else:
    print("general_analysis objects were present before, delete them and start again if needed")

list_contactcutt_frames_conservation=[["45", "50", "4.5ang50pcntframes"],
                                    ["45", "50", "cijcutt0.6_4.5ang50pcntframes"],
                                    ["45", "75", "cijcutt0.6_4.5ang75pcntframes"],
                                      ["6", "75", "6ang75pcntframes"],
                                      ["45", "20:90", "4.5ang20_90pcntframes"],
                                      ["6", "20:90", "6ang20_90pcntframes"]]
jobs = []
desired_numpy_matrices = []
for iterations in list_contactcutt_frames_conservation:
    namesoffile = numpymat+iterations[2]
    # elemenst will be  mainout+"general/contactmaps_numpymatrices/4.5ang50pcntframes"
    desired_numpy_matrices += [namesoffile]
    print (namesoffile)
    if not os.path.isfile(namesoffile+"_resmat"):
        # the contact file doesnt exist, lets create one,
        distcutt = iterations[0] if iterations[0] != '45' else '4.5'
        frames='False' if not replicates else frames
        res=subprocess.check_output(
        ["python", "matrix_prob_counter.py", numpyhas, distcutt, iterations[1], frames, namesoffile])
        print (res)
#sys.exit()
jobs=[]
for numpycmap in desired_numpy_matrices:
    numpycmapob=numpycmap+"_resmat"+"_CMAP.ob"
    if not os.path.isfile(numpycmapob):
        p=multiprocessing.Process(target = numpycmaprunner, args = (
            ultrsmalldcd, cijob, numpycmap+"_resmat", numpycmap+"_resmat"))
        p.daemon=True
        jobs.append(p)
        p.start()
if jobs:
    for p in jobs:
        p.join()
# fine and working
if "rep" in systemtype:
    print("ending replicates before further mess")
    sys.exit()
jobs=[]
analyslis_jobs=[]
for ind, val in enumerate(list_contactcutt_frames_conservation):
    netobfileraw=netobout+val[2]+"net_max_mod.ob"
    netobfilemodref=netobout+val[2]+"net_opt_mod.ob"
    cmapob=desired_numpy_matrices[ind]+"_resmat"+"_CMAP.ob"
    netobfilerawamend=netobout+"amendment"+val[2]+"net_max_mod.ob"
    netobfilemodrefamend=netobout+"amendment"+val[2]+"net_opt_mod.ob"

    netobfilemodref_MEAN=netobout+val[2]+"net_opt_mod_MEAN.ob"
    netobfilemodrefamend_MEAN=netobout+"amendment"+val[2]+"net_opt_mod_MEAN.ob"
    
    if not (os.path.isfile(netobfileraw) and os.path.isfile(netobfilemodref)):
        p=multiprocessing.Process(target = netobcreator, args = (
            cijob, cmapob, val[0], netobout+val[2]))
        p.daemon=True
        jobs.append(p)
        p.start()
        analyslis_jobs+=[["python","reprsentingnetob_31mar2020.py", netobfileraw, "cij", "./", "./", "1"]]
        analyslis_jobs+=[["python","reprsentingnetob_31mar2020.py", netobfilemodref, "cij", "./", "./", "1"]]
        analyslis_jobs+=[["python","reprsentingnetob_31mar2020.py", netobfilerawamend, "cij", "./", "./", "0"]]
        analyslis_jobs+=[["python","reprsentingnetob_31mar2020.py", netobfilemodrefamend, "cij", "./", "./", "0"]]
        
        analyslis_jobs+=[["python","reprsentingnetob_31mar2020.py", netobfilemodref_MEAN, "cij", "./", "./", "1"]]
        analyslis_jobs+=[["python","reprsentingnetob_31mar2020.py", netobfilemodrefamend_MEAN, "cij", "./", "./", "0"]]

if jobs:
    for p in jobs:
        p.join()

if analyslis_jobs:
    for job in analyslis_jobs:
        subprocess.check_output(job)
# represnettaing the objects now

#python comprehensive_python_handler_CIJ.py ensemble ../../source/ensy321a/productiontraj.dcd ../../source/ensy321a/renumber_raw.pdb ../../derived/ensy321a/ ../../source/ensy321a/ultrasmall.dcd ../../source/ensy321a/contacts_closest_heavy_nonlocal.hash 0 20
#python comprehensive_python_handler_CIJ.py ensemble ../../source/enswt/productiontraj.dcd ../../source/enswt/renumber_raw.pdb ../../derived/enswt/ ../../source/enswt/ultrasmall.dcd ../../source/enswt/contacts_closest_heavy_nonlocal.hash 0 20
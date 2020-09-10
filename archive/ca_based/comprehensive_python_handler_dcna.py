import subprocess
import os
from shutil import rmtree
import sys
if len(sys.argv) != 8:
    print(len(sys.argv))
    print("""
    Please enter the correct command line arguements:
    1:ultrasmalldcd 2:pdb 3:outputdir(mutual of wildtyep and mutant)
    4: mutant resmat 5:wildtyperesmat 6 WT prob 7 MT probfile
    """)
    sys.exit()


def createdir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


prog, ultrsmalldcd, pdb, mainout, mat1, mat2, probwt, probmt = sys.argv

if not os.path.isdir(mainout):
    createdir(mainout)


consensusname = mainout+"consensus_matrix"
createsonsensus = subprocess.check_output(
    ["python", "consensus_matrix_of_two.py", mat1, mat2, consensusname])

resmatname = consensusname+"_resmat"

drawdcna = subprocess.check_output(
    ["Rscript", "dcnaapproach.r", ultrsmalldcd, pdb, mainout, resmatname])
networkopt = mainout+'net_opt_mod.ob'
networkmax = mainout+'net_max_mod.ob'

networkoptamend = mainout+'amendment'+'net_opt_mod.ob'
networkmaxamend = mainout+'amendment'+'net_max_mod.ob'
repmax = subprocess.check_output(
    ["python", "reprsentingnetob_31mar2020.py", networkmax, "dcna", probwt, probmt, "1"])
repmaxamend = subprocess.check_output(
    ["python", "reprsentingnetob_31mar2020.py", networkmaxamend, "dcna", probwt, probmt, "0"])

repopt = subprocess.check_output(
    ["python", "reprsentingnetob_31mar2020.py", networkopt, "dcna", probwt, probmt, "1"])
repoptamend = subprocess.check_output(
    ["python", "reprsentingnetob_31mar2020.py", networkoptamend, "dcna", probwt, probmt, "0"])

# above was for nonlocal closesst heavy,,, doing below for cloest heavy  all
# python comprehensive_python_handler_dcna.py ../../source/enswt/ultrasmall.dcd ../../source/enswt/renumber_raw.pdb  ../../derived/ensembledcna/4.5ang90pcnt_closheavy_allcontacts/ ../../derived/ensy321a/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes_.npyob.npy ../../derived/enswt/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang90pcntframes_.npyob.npy ../../derived/enswt/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang_prob ../../derived/ensy321a/general/contactmaps_numpymatrices/allcontacts_closest_heavy_all_contacts_4.5ang_prob

# dcna original approach
# maek consensus matrix of three and consenssu prob before going below, (need microsecond reuslt for good results)
# python comprehensive_python_handler_dcna.py ../../source/enswt/ultrasmall.dcd ../../source/enswt/renumber_raw.pdb ../../derived/ensembledcna/dcna20_90_strict/networks/4.5ang90pcnt/ ../../derived/ensembledcna/dcna20_90_strict/y321a_4.5ang90consensus_.npyob.npy ../../derived/ensembledcna/dcna20_90_strict/WT_4.5ang90consensus_.npyob.npy ../../derived/ensembledcna/dcna20_90_strict/WT_4.5ang_consensus ../../derived/ensembledcna/dcna20_90_strict/y321a_4.5ang_consensus
# python comprehensive_python_handler_dcna.py ../../source/enswt/ultrasmall.dcd ../../source/enswt/renumber_raw.pdb ../../derived/ensembledcna/dcna20_90_strict/networks/6ang90pcnt/ ../../derived/ensembledcna/dcna20_90_strict/y321a_6ang90consensus_.npyob.npy ../../derived/ensembledcna/dcna20_90_strict/WT_6ang90consensus_.npyob.npy ../../derived/ensembledcna/dcna20_90_strict/WT_6ang_consensus ../../derived/ensembledcna/dcna20_90_strict/y321a_6ang_consensus

# dcna runs 27 may
#

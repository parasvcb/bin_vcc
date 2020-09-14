import os
import subprocess
import numpy as np
import sys
import re
from pathlib import Path
import module_contacts
if len(sys.argv) != 7:
    print("Please enter correct cmd souercedd rep1add rep2add rep3add combinedtraj cijcutt processingdir")
    sys.exit()
prog, add1, add2, add3, concattraj, cuttcij, outloc = sys.argv
# now source add will have source diles,a dn last appebd of them (or name), will be used to create directory in
# processingdir and files like rmsf rmsd will be serahed over there

# python process.py ../../ready_analysis/Y321F_r1/ ../../ready_analysis/Y321F_r2/ ../../ready_analysis/Y321F_r3/ ../../ready_analysis/combined/ 0.6 ../../output_process/


def checkandcreatedir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


pdbtype = 'coordmini.pdb'
pdbtype = 'renumber_raw.pdb'


rep1name = os.path.basename(os.path.normpath(add1))
rep2name = os.path.basename(os.path.normpath(add2))
rep3name = os.path.basename(os.path.normpath(add3))

pro_rep1 = os.path.join(outloc, rep1name)
pro_rep2 = os.path.join(outloc, rep2name)
pro_rep3 = os.path.join(outloc, rep3name)

proreplis = {pro_rep1: add1, pro_rep2: add2, pro_rep3: add3}
checkandcreatedir(pro_rep1)
checkandcreatedir(pro_rep2)
checkandcreatedir(pro_rep3)
# print(proreplis)
# sys.exit()
# contactswill be named as contacts_ca-ca_nonlocal.hash
for i in proreplis:
    parent = proreplis[i]
    contactshas = os.path.join(i, 'contacts_ca-ca_nonlocal.hash')
    contactsmatrix = os.path.join(
        i, 'contacts_ca-ca_nonlocal_matrix_10_75_resmat')
    cijmatrix = os.path.join(i, 'sim_data', 'mymatrixcij_CA.txt')
    checkandcreatedir(os.path.join(i, 'sim_data'))
    dcd = os.path.join(parent, 'productiontraj.dcd')
    pdb = os.path.join(parent, pdbtype)

    if not os.path.isfile(contactshas):
        module_contacts.horse(dcd, pdb, i)
    # print(os.path.isfile(contactshas))
    # print(contactshas)
    if not os.path.isfile(contactsmatrix):
        mat, probhas = module_contacts.matrix_returner_single_cuttoff(
            contactshas, 10, 0.75, upp=False, low=False)
        module_contacts.write_matrix(mat, os.path.join(
            i, 'contacts_ca-ca_nonlocal_matrix_10_75'))
        module_contacts.writeprobfile(probhas, os.path.join(
            i, 'nonnumbered_contacts_probdiff_10_75.text'))
    if not os.path.isfile(cijmatrix):
        parentpathtemp = os.path.dirname(os.path.normpath(cijmatrix))
        module_contacts.callcij(dcd, pdb, parentpathtemp)

cuttcij = 0.6
consensusFile = os.path.join(outloc, 'CONSENSUS_Cij_MATRIX')
if not os.path.isfile(consensusFile+'_resmat'):

    mat1contact = os.path.join(
        pro_rep1, "contacts_ca-ca_nonlocal_matrix_10_75_.npyob.npy")
    mat2contact = os.path.join(
        pro_rep2, "contacts_ca-ca_nonlocal_matrix_10_75_.npyob.npy")
    mat3contact = os.path.join(
        pro_rep3, "contacts_ca-ca_nonlocal_matrix_10_75_.npyob.npy")

    con1 = np.load(mat1contact)
    con2 = np.load(mat2contact)
    con3 = np.load(mat3contact)

    mat1cij = os.path.join(pro_rep1, "sim_data", "mymatrixcij_CA.txt")
    mat2cij = os.path.join(pro_rep2, "sim_data", "mymatrixcij_CA.txt")
    mat3cij = os.path.join(pro_rep3, "sim_data", "mymatrixcij_CA.txt")

    mat1 = module_contacts.texttomatrix(mat1cij)
    mat2 = module_contacts.texttomatrix(mat2cij)
    mat3 = module_contacts.texttomatrix(mat3cij)

    module_contacts.consensus_matrix(
        mat1, mat2, mat3, con1, con2, con3, cuttcij, consensusFile)

jbcdir = os.path.join(outloc, 'jbcapproach')
checkandcreatedir(jbcdir)
netobdefaultfile = os.path.join(jbcdir, 'net_max_mod.ob')
if not os.path.isfile(netobdefaultfile):
    pdb = os.path.join(add1, pdbtype)
    ultrasmalldcd = os.path.join(add1, 'ultrasmall.dcd')
    consfile = consensusFile+'_resmat'
    module_contacts.createnetob(ultrasmalldcd, pdb, consfile, jbcdir)


# next module will be final and should be ready()
# it will give a preprocessed R verison, that eeds to be restrcutred using python
# each strcutre shluld have community wise data, defaulty edge weight and complete correlations in side
detailedinfofile = os.path.join(jbcdir, 'net_max_moddetailed.text')
cijfilehastext = os.path.join(
    outloc, 'CONSENSUS_Cij_MATRIX_valuesCij.hashtext')
# if not os.path.isfile(detailedinfofile):
if 1:
    tempappend = os.path.join(jbcdir, 'net_max_mod')
    rawFile = tempappend+"rawfile_cijvalues_communitywisefromR.text"
    anycijob = os.path.join(pro_rep3, "sim_data", "CIJ_CA.ob")
    resmatconsensus = consensusFile+'_resmat'
    qualifierTyperawFile = consensusFile+'_qualifierType.hashtext'
    hasmem = module_contacts.get_members(
        netobdefaultfile, tempappend)
    rawfileData = module_contacts.quickRparse(
        tempappend+'.ob', anycijob, resmatconsensus,  rawFile)
    # print('yes3')
    #print(cijfilehastext, 'here')
    # print(hasmem)
    module_contacts.processfinaldata(
        rawfileData, cijfilehastext, hasmem, qualifierTyperawFile, detailedinfofile)

if 0:
    # for i in proreplis:
    #     parent = proreplis[i]
    #     ultrasmalldcd = os.path.join(parent, 'ultrasmall.dcd')
    #     dcd = os.path.join(parent, 'productiontraj.dcd')
    #     pdb = os.path.join(parent, pdbtype)
    #     module_contacts.generalAnalysis(dcd, pdb, i)
    conname = os.path.basename(os.path.normpath(concattraj))
    resconname = os.path.join(outloc, conname)
    checkandcreatedir(resconname)
    ultrasmalldcd = os.path.join(concattraj, 'ultrasmall.dcd')
    dcd = os.path.join(concattraj, 'productiontraj.dcd')
    pdb = os.path.join(concattraj, pdbtype)
    module_contacts.generalAnalysis(dcd, pdb, resconname)
print('done')


# python process.py ../../../WT_series/ready_analysis/WT_r1/ ../../../WT_series/ready_analysis/WT_r2/ ../../../WT_series/ready_analysis/WT_r3/ ../../../WT_series/ready_analysis/combined/ 0.6 ../../../WT_series/output_process/
# python process.py ../../../Y321A_series/ready_analysis/Y321A_r4as1/ ../../../Y321A_series/ready_analysis/Y321A_r2/ ../../../Y321A_series/ready_analysis/Y321A_r3/ ../../../Y321A_series/ready_analysis/combined/ 0.6 ../../../Y321A_series/output_process/
# python process.py ../../ready_analysis/Y321F_r1/ ../../ready_analysis/Y321F_r2/ ../../ready_analysis/Y321F_r3/ ../../ready_analysis/combined/ 0.6 ../../output_process/

# python process.py ../../../WT_series/ready_analysis/WT_r1/ ../../../WT_series/ready_analysis/WT_r2/ ../../../WT_series/ready_analysis/WT_r3/ ../../../WT_series/ready_analysis/combined/ 0.6 ../../../WT_series/output_process_renumberpdb_aligned/
# python process.py ../../ready_analysis/Y321F_r1/ ../../ready_analysis/Y321F_r2/ ../../ready_analysis/Y321F_r3/ ../../ready_analysis/combined/ 0.6 ../../output_process_renumberpdb_aligned/
# python process.py ../../../Y321A_series/ready_analysis/Y321A_r4as1/ ../../../Y321A_series/ready_analysis/Y321A_r2/ ../../../Y321A_series/ready_analysis/Y321A_r3/ ../../../Y321A_series/ready_analysis/combined/ 0.6 ../../../Y321A_series/output_process_renumberpdb_aligned/

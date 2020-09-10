import os
import subprocess
import numpy as np
import sys
import re
from pathlib import Path
import module_contacts
if len(sys.argv) != 7:
    print("Please enter correct cmd souercedd rep1add rep2add rep3add cijcutt processingdir")
    sys.exit()
prog, add1, add2, add3, concattraj, cuttcij, outloc = sys.argv
# now source add will have source diles,a dn last appebd of them (or name), will be used to create directory in
# processingdir and files like rmsf rmsd will be serahed over there


def checkandcreatedir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


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

# contactswill be named as contacts_ca-ca_nonlocal.hash
for i in proreplis:
    parent = proreplis[i]
    contactshas = os.path.join(i, 'contacts_ca-ca_nonlocal.hash')
    contactsmatrix = os.path.join(
        i, 'contacts_ca-ca_nonlocal_matrix_10_75_resmat')
    cijmatrix = os.path.join(i, 'sim_data', 'mymatrixcij_CA.txt')
    dcd = os.path.join(parent, 'productiontraj.dcd')
    pdb = os.path.join(parent, 'coordmini.pdb')

    if not os.path.isfile(contactshas):
        module_contacts.horse(dcd, pdb, i)
    if not os.path.isfile(contactsmatrix):
        mat, probhas = module_contacts.matrix_returner_single_cuttoff(contactshas, 10, 0.75, upp=False, low=False):
        module_contacts.write_matrix(mat, os.path.join(
            i, 'contacts_ca-ca_nonlocal_matrix_10_75'))
        module_contacts.writeprobfile(probhas, os.path.join(
            i, 'nonnumbered_contacts_probdiff_10_75.text'))
    if not os.path.isfile(cijmatrix):
        module.callcij(dcd, pdb, Path(cijmatrix).parent.name)

prog, add1, add2, add3, cuttcij, postappend, resfile = sys.argv
cuttcij = 0.6
consensusFile = os.path.join(outloc, 'CONSENSUS_Cij_MATRIX')
if not os.path.isfile(consensusFile+'_resmat')):

    mat1contact=os.path.join(
        pro_rep1, "contacts_ca-ca_nonlocal_matrix_10_75_.npyob.npy")
    mat2contact=os.path.join(
        pro_rep2, "contacts_ca-ca_nonlocal_matrix_10_75_.npyob.npy")
    mat3contact=os.path.join(
    pro_rep3, "contacts_ca-ca_nonlocal_matrix_10_75_.npyob.npy")

    con1 = np.load(mat1contact)
    con2 = np.load(mat2contact)
    con3 = np.load(mat3contact)

    mat1cij = os.path.join(pro_rep1, "sim_data/mymatrixcij_CA.txt")
    mat2cij = os.path.join(pro_rep2, "sim_data/mymatrixcij_CA.txt")
    mat3cij = os.path.join(pro_rep3, "sim_data/mymatrixcij_CA.txt")

    mat1 = module_contacts.texttomatrix(mat1cij)
    mat2 = module_contacts.texttomatrix(mat2cij)
    mat3 = module_contacts.texttomatrix(mat3cij)

    module_contacts.consensus_matrix(mat1, mat2, mat3, con1, con2, con3, cuttcij, consensusFile)

jbcdir=os.path.join(outloc,'jbcapproach')
checkandcreatedir(jbcdir)
netobdefaultfile=os.path.join(jbcdir,'net_max_mod.ob')
if not os.path.isfile(netobdefaultfile):
    pdb=os.path.join(proreplis[1], 'coordmini.pdb')
    ultrasmalldcd=os.path.join(proreplis[1], 'ultrasmall.dcd')
    consfile=os.path.isfile(consensusFile+'_resmat')
    module_contacts.createnetob(ultrasmalldcd,pdb,consfile,jbcdir)


# next module will be final and should be ready()
# it will give a preprocessed R verison, that eeds to be restrcutred using python
# each strcutre shluld have community wise data, defaulty edge weight and complete correlations in side
detailedinfofile=os.path.join(jbcdir,'net_max_moddetailed.text')
cijfilehastext=os.path.join(outloc, 'CONSENSUS_Cij_MATRIX_valuesCij.hashtext')
if not os.path.isfile(detailedinfofile):
    tempappend=os.path.join(jbcdir+'net_max_mod')
    rawFile=tempappend+"rawfile_cijvalues_communitywisefromR.text"
    hasmem=module_contacts.get_members(netobdefaultfile,tempappend)
    rawfileData=module_contacts.postprocessing(tempappend+'.ob',rawFile)
    module_contacts.processfinaldata(rawfileData,cijfilehastext,hasmem)

if general_analysis:
    for i in proreplis:
        parent = proreplis[i]
        ultrasmalldcd = os.path.join(parent, 'ultrasmall.dcd') 
        dcd = os.path.join(parent, 'productiontraj.dcd')
        pdb = os.path.join(parent, 'coordmini.pdb')
        module_contacts.generalAnalysis(dcd,pdb,i)
    conname = os.path.basename(os.path.normpath(concattraj))
    resconname = os.path.join(outloc, conname)
    checkandcreatedir(resconname)
    ultrasmalldcd = os.path.join(concattraj, 'ultrasmall.dcd') 
    dcd = os.path.join(concattraj, 'productiontraj.dcd')
    pdb = os.path.join(concattraj, 'coordmini.pdb')
    module_contacts.generalAnalysis(dcd,pdb,resconname)
print('done')



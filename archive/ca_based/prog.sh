python comprehensive_python_handler_CIJ.py ensemble ../../source/enswt/productiontraj.dcd ../../source/enswt/renumber_raw.pdb ../../derived/enswt/ ../../source/enswt/ultrasmall.dcd ../../source/enswt/contacts.hash 0 2
python comprehensive_python_handler_CIJ.py ensemble ../../source/ensy321a/productiontraj.dcd ../../source/ensy321a/renumber_raw.pdb ../../derived/ensy321a/ ../../source/ensy321a/ultrasmall.dcd ../../source/ensy321a/contacts.hash 0 2

# python comprehensive_python_handler_CIJ.py replicates_0:10000 ../../source/wt_r1/productiontraj.dcd ../../source/wt_r1/renumber_raw.pdb ../../derived/wt_r1/ ../../source/wt_r1/ultrasmall.dcd ../../source/enswt/contacts.hash 0 2
# python comprehensive_python_handler_CIJ.py replicates_10000:20000 ../../source/wt_r2/productiontraj.dcd ../../source/wt_r2/renumber_raw.pdb ../../derived/wt_r2/ ../../source/wt_r2/ultrasmall.dcd ../../source/enswt/contacts.hash 0 2
# python comprehensive_python_handler_CIJ.py replicates_20000:30000 ../../source/wt_r3/productiontraj.dcd ../../source/wt_r3/renumber_raw.pdb ../../derived/wt_r3/ ../../source/wt_r3/ultrasmall.dcd ../../source/enswt/contacts.hash 0 2

# python comprehensive_python_handler_CIJ.py replicates_0:10000 ../../source/y321a_r4/productiontraj.dcd ../../source/y321a_r4/renumber_raw.pdb ../../derived/y321a_r4/ ../../source/y321a_r4/ultrasmall.dcd ../../source/ensy321a/contacts.hash 0 2
# python comprehensive_python_handler_CIJ.py replicates_10000:20000 ../../source/y321a_r2/productiontraj.dcd ../../source/y321a_r2/renumber_raw.pdb ../../derived/y321a_r2/ ../../source/y321a_r2/ultrasmall.dcd ../../source/ensy321a/contacts.hash 0 2
# python comprehensive_python_handler_CIJ.py replicates_20000:30000 ../../source/y321a_r3/productiontraj.dcd ../../source/y321a_r3/renumber_raw.pdb ../../derived/y321a_r3/ ../../source/y321a_r3/ultrasmall.dcd ../../source/ensy321a/contacts.hash 0 2




# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 4.5 90 False derived/enswt/general/contactmaps_numpymatrices/4.5ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 4.5 90  False derived/ensy321a/general/contactmaps_numpymatrices/4.5ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 4.5 90 0:10000  derived/wt_r1/general/contactmaps_numpymatrices/4.5ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 4.5 90 10000:20000 derived/wt_r2/general/contactmaps_numpymatrices/4.5ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 4.5 90  20000:30000 derived/wt_r3/general/contactmaps_numpymatrices/4.5ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 4.5 90  0:10000    derived/y321a_r4/general/contactmaps_numpymatrices/4.5ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 4.5 90  10000:20000 derived/y321a_r2/general/contactmaps_numpymatrices/4.5ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 4.5 90 20000:30000 derived/y321a_r3/general/contactmaps_numpymatrices/4.5ang90pcntframes 1



# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 6 90 False derived/enswt/general/contactmaps_numpymatrices/6ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 6 90  False derived/ensy321a/general/contactmaps_numpymatrices/6ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 6 90 0:10000  derived/wt_r1/general/contactmaps_numpymatrices/6ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 6 90 10000:20000 derived/wt_r2/general/contactmaps_numpymatrices/6ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 6 90  20000:30000 derived/wt_r3/general/contactmaps_numpymatrices/6ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 6 90  0:10000    derived/y321a_r4/general/contactmaps_numpymatrices/6ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 6 90  10000:20000 derived/y321a_r2/general/contactmaps_numpymatrices/6ang90pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 6 90 20000:30000 derived/y321a_r3/general/contactmaps_numpymatrices/6ang90pcntframes 1


# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 10 75 False derived/enswt/general/contactmaps_numpymatrices/10ang75pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 10 75  False derived/ensy321a/general/contactmaps_numpymatrices/10ang75pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 10 75 0:10000  derived/wt_r1/general/contactmaps_numpymatrices/10ang75pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 10 75 10000:20000 derived/wt_r2/general/contactmaps_numpymatrices/10ang75pcntframes 1
# python bin/ca_based/prob_counter.py source/enswt/contacts.hash 10 75  20000:30000 derived/wt_r3/general/contactmaps_numpymatrices/10ang75pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 10 75  0:10000    derived/y321a_r4/general/contactmaps_numpymatrices/10ang75pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 10 75  10000:20000 derived/y321a_r2/general/contactmaps_numpymatrices/10ang75pcntframes 1
# python bin/ca_based/prob_counter.py source/ensy321a/contacts.hash 10 75 20000:30000 derived/y321a_r3/general/contactmaps_numpymatrices/10ang75pcntframes 1


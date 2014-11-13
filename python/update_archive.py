import h5py
import sys
import numpy
import subprocess

if len(sys.argv) < 2:
  print "Usage: python update_archive.py old_archive"
  sys.exit()

print """
This script is an attempt to update your archive to TRIQS 1.2.
Please keep a copy of your old archive as this script is
** not guaranteed ** to work for your archive.
If you encounter any problem please report it on github!
"""

filename = sys.argv[1]
A = h5py.File(filename)

old_to_new = {'SumK_LDA':'lda_input', 'SumK_LDA_ParProj':'lda_parproj_input', 
 'SymmCorr':'lda_symmcorr_input', 'SymmPar':'lda_symmpar_input', 'SumK_LDA_Bands':'lda_bands_input'}

for old, new in old_to_new.iteritems():
    if old not in A.keys(): continue
    print "Changing %s to %s ..."%(old, new)
    A.copy(old,new)
    del(A[old])

move_to_output = ['gf_struct_solver','map_inv','map',
                  'chemical_potential','dc_imp','dc_energ','deg_shells',
                  'h_field']
for obj in move_to_output:
    if obj in A['lda_input'].keys():
       if not 'lda_output' in A: A.create_group('lda_output')
       print "Moving %s to lda_output ..."%obj
       A.copy('lda_input/'+obj,'lda_output/'+obj)
       del(A['lda_input'][obj])

A.close()

# Repack to reclaim disk space
retcode = subprocess.call(["h5repack","-i%s"%filename, "-otemphgfrt.h5"])
if retcode != 0:
    print "h5repack failed!"
else:
    subprocess.call(["mv","-f","temphgfrt.h5","%s"%filename])

from pytriqs.archive import HDFArchive
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

def det_shell_equivalence(lst):
    corr_to_inequiv = [0 for i in range(len(lst))]
    inequiv_to_corr = [0]
    n_inequiv_shells = 1
    tmp = [ lst[0][1:3] ]
    if (len(lst)>1):
        for i in range(len(lst)-1):
            fnd = False
            for j in range(n_inequiv_shells):
                if (tmp[j]==lst[i+1][1:3]):
                    fnd = True
                    corr_to_inequiv[i+1] = j
            if (fnd==False):
                corr_to_inequiv[i+1] = n_inequiv_shells
                n_inequiv_shells += 1
                tmp.append( lst[i+1][1:3] )
                inequiv_to_corr.append(i+1)
    return [n_inequiv_shells, corr_to_inequiv, inequiv_to_corr]

### Main ###

filename = sys.argv[1]
A = h5py.File(filename)

# Rename groups
old_to_new = {'SumK_LDA':'lda_input', 'SumK_LDA_ParProj':'lda_parproj_input', 
 'SymmCorr':'lda_symmcorr_input', 'SymmPar':'lda_symmpar_input', 'SumK_LDA_Bands':'lda_bands_input'}

for old, new in old_to_new.iteritems():
    if old not in A.keys(): continue
    print "Changing %s to %s ..."%(old, new)
    A.copy(old,new)
    del(A[old])

# Move output items from lda_input to lda_output
move_to_output = ['gf_struct_solver','map_inv','map',
                  'chemical_potential','dc_imp','dc_energ','deg_shells',
                  'h_field']
for obj in move_to_output:
    if obj in A['lda_input'].keys():
       if not 'lda_output' in A: A.create_group('lda_output')
       print "Moving %s to lda_output ..."%obj
       A.copy('lda_input/'+obj,'lda_output/'+obj)
       del(A['lda_input'][obj])

# Add shell equivalency quantities
B = A['lda_input']
corr_shells = HDFArchive(filename,'r')['lda_input']['corr_shells']
equiv_shell_info = det_shell_equivalence(corr_shells)
B['n_inequiv_shells'] = equiv_shell_info[0]
B['corr_to_inequiv'] = equiv_shell_info[1]
B['inequiv_to_corr'] = equiv_shell_info[2]

# Rename variables
groups = ['lda_symmcorr_input','lda_symmpar_input']

for group in groups:
    if group not in A.keys(): continue
    print "Changing n_s to n_symm ..."
    A[group].move('n_s','n_symm')

A.close()

# Repack to reclaim disk space
retcode = subprocess.call(["h5repack","-i%s"%filename, "-otemphgfrt.h5"])
if retcode != 0:
    print "h5repack failed!"
else:
    subprocess.call(["mv","-f","temphgfrt.h5","%s"%filename])

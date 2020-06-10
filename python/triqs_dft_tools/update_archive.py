from h5 import HDFArchive
import h5py
import sys
import numpy
import subprocess

if len(sys.argv) < 2:
    print("Usage: python update_archive.py old_archive [v1.0|v1.2]")
    sys.exit()

print("""
This script is an attempt to update your archive to TRIQS 1.2.
Please keep a copy of your old archive as this script is
** not guaranteed ** to work for your archive.
If you encounter any problem please report it on github!
""")


def convert_shells(shells):
    shell_entries = ['atom', 'sort', 'l', 'dim']
    return [{name: int(val) for name, val in zip(shell_entries, shells[ish])} for ish in range(len(shells))]


def convert_corr_shells(corr_shells):
    corr_shell_entries = ['atom', 'sort', 'l', 'dim', 'SO', 'irep']
    return [{name: int(val) for name, val in zip(corr_shell_entries, corr_shells[icrsh])} for icrsh in range(len(corr_shells))]


def det_shell_equivalence(corr_shells):
    corr_to_inequiv = [0 for i in range(len(corr_shells))]
    inequiv_to_corr = [0]
    n_inequiv_shells = 1

    if len(corr_shells) > 1:
        inequiv_sort = [corr_shells[0]['sort']]
        inequiv_l = [corr_shells[0]['l']]
        for i in range(len(corr_shells) - 1):
            is_equiv = False
            for j in range(n_inequiv_shells):
                if (inequiv_sort[j] == corr_shells[i + 1]['sort']) and (inequiv_l[j] == corr_shells[i + 1]['l']):
                    is_equiv = True
                    corr_to_inequiv[i + 1] = j
            if is_equiv == False:
                corr_to_inequiv[i + 1] = n_inequiv_shells
                n_inequiv_shells += 1
                inequiv_sort.append(corr_shells[i + 1]['sort'])
                inequiv_l.append(corr_shells[i + 1]['l'])
                inequiv_to_corr.append(i + 1)

    return n_inequiv_shells, corr_to_inequiv, inequiv_to_corr


### Main ###

filename = sys.argv[1]
if len(sys.argv) > 2:
    from_v = sys.argv[2]
else:  # Assume updating an old v1.0 script
    from_v = 'v1.0'
A = h5py.File(filename)

# Rename groups
old_to_new = {'SumK_LDA': 'dft_input', 'SumK_LDA_ParProj': 'dft_parproj_input',
              'SymmCorr': 'dft_symmcorr_input', 'SymmPar': 'dft_symmpar_input', 'SumK_LDA_Bands': 'dft_bands_input'}

for old, new in old_to_new.items():
    if old not in list(A.keys()):
        continue
    print("Changing %s to %s ..." % (old, new))
    A.copy(old, new)
    del(A[old])

# Move output items from dft_input to user_data
move_to_output = ['chemical_potential', 'dc_imp', 'dc_energ']
for obj in move_to_output:
    if obj in list(A['dft_input'].keys()):
        if 'user_data' not in A:
            A.create_group('user_data')
        print("Moving %s to user_data ..." % obj)
        A.copy('dft_input/' + obj, 'user_data/' + obj)
        del(A['dft_input'][obj])
# Delete obsolete quantities
to_delete = ['gf_struct_solver', 'map_inv', 'map', 'deg_shells', 'h_field']
for obj in to_delete:
    if obj in list(A['dft_input'].keys()):
        del(A['dft_input'][obj])

if from_v == 'v1.0':
    # Update shells and corr_shells to list of dicts
    shells_old = HDFArchive(filename, 'r')['dft_input']['shells']
    corr_shells_old = HDFArchive(filename, 'r')['dft_input']['corr_shells']
    shells = convert_shells(shells_old)
    corr_shells = convert_corr_shells(corr_shells_old)
    del(A['dft_input']['shells'])
    del(A['dft_input']['corr_shells'])
    A.close()
    # Need to use HDFArchive for the following
    HDFArchive(filename, 'a')['dft_input']['shells'] = shells
    HDFArchive(filename, 'a')['dft_input']['corr_shells'] = corr_shells
    A = h5py.File(filename)

# Add shell equivalency quantities
if 'n_inequiv_shells' not in A['dft_input']:
    equiv_shell_info = det_shell_equivalence(corr_shells)
    A['dft_input']['n_inequiv_shells'] = equiv_shell_info[0]
    A['dft_input']['corr_to_inequiv'] = equiv_shell_info[1]
    A['dft_input']['inequiv_to_corr'] = equiv_shell_info[2]

# Rename variables
groups = ['dft_symmcorr_input', 'dft_symmpar_input']
for group in groups:
    if group not in list(A.keys()):
        continue
    if 'n_s' not in A[group]:
        continue
    print("Changing n_s to n_symm ...")
    A[group].move('n_s', 'n_symm')
    # Convert orbits to list of dicts
    orbits_old = HDFArchive(filename, 'r')[group]['orbits']
    orbits = convert_corr_shells(orbits_old)
    del(A[group]['orbits'])
    A.close()
    HDFArchive(filename, 'a')[group]['orbits'] = orbits
    A = h5py.File(filename)

groups = ['dft_parproj_input']
for group in groups:
    if group not in list(A.keys()):
        continue
    if 'proj_mat_pc' not in A[group]:
        continue
    print("Changing proj_mat_pc to proj_mat_all ...")
    A[group].move('proj_mat_pc', 'proj_mat_all')

A.close()

# Repack to reclaim disk space
retcode = subprocess.call(["h5repack", "-i%s" % filename, "-otemphgfrt.h5"])
if retcode != 0:
    print("h5repack failed!")
else:
    subprocess.call(["mv", "-f", "temphgfrt.h5", "%s" % filename])

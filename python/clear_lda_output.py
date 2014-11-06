import h5py
import sys
import subprocess

if len(sys.argv) < 2:
  print "Usage: python clear_lda_output.py archive"
  sys.exit()

print """
This script is to remove any SumkLDA generated output from the h5 archive
and to restore it to the original post-converter state.
"""

filename = sys.argv[1]
A = h5py.File(filename)
del(A["lda_output"])
A.close()

# Repack to reclaim disk space
retcode = subprocess.call(["h5repack","-i%s"%filename, "-otemphgfrt.h5"])
if retcode != 0:
    print "h5repack failed!"
else:
    subprocess.call(["mv","-f","temphgfrt.h5","%s"%filename])

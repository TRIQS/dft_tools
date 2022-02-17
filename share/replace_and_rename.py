#!/usr/bin/env python2

import sys
import os
import glob

if len(sys.argv) != 2:
    print("Please pass the application name")
    sys.exit()

app_name = str(sys.argv[1]).lower()
capital_name = app_name.upper()

# Move app4triqs directories if necessary
if os.path.isdir("c++/app4triqs"): os.rename("c++/app4triqs", "c++/" + app_name)
if os.path.isdir("python/app4triqs"): os.rename("python/app4triqs", "python/" + app_name)

# Ignore these files
ignore_lst = [".git/", "replace_and_rename.py", "squash_history.sh"]

# Find the root directory of app4triqs
app4triqs_root = os.path.abspath(os.path.dirname(__file__) + "/..")

# Recurse over all subdirectories and files
for root, dirs, files in os.walk(app4triqs_root):

    for fname in files:
        fpath = os.path.join(root, fname)

        # Ignore certain files / directories
        if any(it in fpath for it in ignore_lst): continue

        if os.path.isfile(fpath):
            # Rename files containing app4triqs in their filename
            if "app4triqs" in fname:
                new_fpath = os.path.join(root, fname.replace("app4triqs", app_name))
                os.rename(fpath, new_fpath)
                fpath = new_fpath

            # Replace app4triqs and APP4TRIQS in all files
            with open(fpath, 'r') as f:
                s = f.read()
            if "app4triqs" in s or "APP4TRIQS" in s:
                with open(fpath, 'w') as f:
                    f.write(s.replace("app4triqs", app_name).replace("APP4TRIQS", capital_name))

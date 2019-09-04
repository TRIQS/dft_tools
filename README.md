[![Build Status](https://travis-ci.org/TRIQS/app4triqs.svg?branch=unstable)](https://travis-ci.org/TRIQS/app4triqs)

# app4triqs - A skeleton for a TRIQS application

Initial Setup
-------------

**Caution**: The following instructions require the `util-linux` rename command.
Please confirm that you have the right version by running `rename --version`.
For the Perl rename command see instructions in the following section.

To adapt this skeleton for a new TRIQS application, the following steps are necessary:

* Create a repository, e.g. https://github.com/myuser/mynewapp

* Run the following commands in order after replacing myuser, mynewapp and MYNEWAPP accordingly

```bash
github_username=myuser
app_name=mynewapp
capital_name=MYNEWAPP

git clone https://github.com/triqs/app4triqs --branch unstable ${app_name}
cd ${app_name}
git reset $(git commit-tree HEAD\^{tree} -m "Create ${app_name} from github.com/triqs/app4triqs skeleton")
git merge --allow-unrelated-histories -s ours HEAD@{1} -m "Track app4triqs skeleton"
find . -type f | grep -v .git | xargs sed -i 's/app4triqs/${app_name}/g; s/APP4TRIQS/${capital_name}/g'
find . -type d | grep -v .git | xargs rename app4triqs ${app_name}
find . -type f | grep -v .git | xargs rename app4triqs ${app_name}
git add -A && git commit -m "Adjust app4triqs skeleton for ${app_name}"
git remote set-url origin https://github.com/${github_username}/${app_name}
git remote add app4triqs_remote https://github.com/triqs/app4triqs
git remote update && git remote prune origin
```

You can now push to your github repository

```bash
git push origin unstable
```

### Perl rename command ###

If you are using the Perl-based rename command you will need to 

```bash
find . -type d | grep -v .git | xargs rename 's/app4triqs/${app_name}/'
find . -type f | grep -v .git | xargs rename 's/app4triqs/${app_name}/'
```

### Github SSH interface ###

If you prefer to use the SSH interface to the remote repository,
replace the http link accordingly

```
https://github.com/myuser/mynewapp   -->   git@github.com:myuser/mynewapp
```

### Merging app4triqs skeleton updates ###

You can merge future changes to app4triqs into your project with the following commands

```bash
git remote update
git merge app4triqs_remote -m "Merge latest app4triqs skeleton changes"
```

If you should encounter any conflicts resolve them and `git commit`.

Getting Started
---------------

After setting up your application as described above you should customize the following files and directories
according to your needs (replace app4triqs in the following by the name of your application)

* In the `c++/app4triqs` subdirectory adjust the example files `app4triqs.hpp` and `app4triqs.cpp` or add your own source files.
* In the `test/c++` subdirectory adjust the example test `basic.cpp` or add your own tests.
* In the `python/app4triqs` subdirectory add your Python source files.
  Be sure to remove the `app4triqs_module_desc.py` file unless you want to generate a Python module from your C++ source code.
* In the `test/c++` subdirectory adjust the example test `basic.cpp` or add your own tests.
* The build and install process is identical to the one outline [here](https://triqs.github.io/app4triqs/unstable/install.html).

### Optional ###
----------------
* If you want to wrap C++ classes and/or functions provided in the `c++/app4triqs/app4triqs.hpp` rerun the `c++2py` tool with
```bash
c++2py -r app4triqs_module_desc.py
```
* Add your email address to the bottom section of `Jenkinsfile` for Jenkins CI notification emails
```
End of build log:
\${BUILD_LOG,maxLines=60}
    """,
    to: 'user@domain.org',
```

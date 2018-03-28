# app4triqs

An example application using cpp2py and triqs
---------------------------------------------

To use this skeleton for a new triqs application, the following steps are necessary:

* Create a repository, e.g. https://github.com/myuser/mynewapp

* Run the following commands (replacing myuser and mynewapp accordingly)

```bash
git clone https://github.com/triqs/app4triqs --branch unstable mynewapp
cd mynewapp
git remote set-url https://github.com/myuser/mynewapp
find . -type f | grep -v .git | xargs sed -i 's/app4triqs/mynewapp/g; s/APP4TRIQS/MYNEWAPP/g'
find . | grep -v .git | xargs rename app4triqs mytriqsapp
git add -A
git commit -m "Create mynewapp from github.com/triqs/app4triqs skeleton"
git push
```

If you prefer to use the SSH interface to the remote repository,
replace the http link accordingly

```
https://github.com/myuser/mynewapp   -->   git@github.com:myuser/mynewapp
```

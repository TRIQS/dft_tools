#!/usr/bin/env bash

git reset $(git commit-tree HEAD\^{tree} -m "Initialize project from github.com/triqs/app4triqs@$(git rev-parse --short HEAD)")
git merge --allow-unrelated-histories -s ours HEAD@{1} -m "Track app4triqs skeleton"
git remote rm origin
git remote add app4triqs_remote https://github.com/triqs/app4triqs

#!/usr/bin/env bash

if [ $# -ne 1 ]; then
  echo "Please pass the application name"
  exit 1
fi
app_name=$1
capital_name=$(printf '%s' "$1" | awk '{ print toupper($0) }')

# Move app4triqs directories if necessary
[ -d c++/app4triqs ] && mv c++/app4triqs c++/${app_name}
[ -d python/app4triqs ] && mv python/app4triqs python/${app_name}

# Replace app4triqs and APP4TRIQS for our application in all files and filenames
if [ $(uname -s) == Linux ]; then
  find . -type f \
    -not -path "./.git/*" \
    -not -path "*/replace_and_rename.sh" \
    -not -path "*/squash_history.sh" \
    -exec sed -i "s/app4triqs/${app_name}/g; s/APP4TRIQS/${capital_name}/g" {} \;
  find . -type f -not -path "./.git/*" -exec rename app4triqs ${app_name} {} &> /dev/null \;
elif [ $(uname -s) == Darwin ]; then
  LC_CTYPE=C LANG=C find . -type f \
    -not -path "./.git/*" \
    -not -path "*/replace_and_rename.sh" \
    -not -path "*/squash_history.sh" \
    -exec sed -i '' -e "s/app4triqs/${app_name}/g; s/APP4TRIQS/${capital_name}/g" {} \;
  find . -type f -not -path "./.git/*" -exec rename "s/app4triqs/${app_name}/" {} &> /dev/null \;
fi

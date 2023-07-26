#! /bin/bash

git submodule init
git submodule update --remote 

# to update the repo (after changes are made to the parent repo), run:
# > git fetch
# > git pull origin main

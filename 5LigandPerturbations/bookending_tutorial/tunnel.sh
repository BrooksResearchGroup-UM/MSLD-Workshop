#! /usr/bin/bash

# Load modules
module load cgenff
module load mmtsb/mmtsb
module load anaconda/2022
module load openbabel

# Load conda environment with pyCHARMM installed
conda activate pycharmm47

# Make sure you are in a node in the cluster
intra_ip=`ifconfig | grep "inet 192.*" | awk '{print $2;}'`
jupyter-notebook --no-browser --ip $intra_ip --port 49153 --no-script


# On local do something similar to:
# `ssh -N -f -L localhost:9999:gollum152:49153 gollum`
#
# The 9999 is the listening port on remote
#
# The 49153 is the broadcasting port from the cluster
# 
# gollum152 is the node where I am running things from
#
# gollum is the cluster name 

#! /bin/bash

sbatch --time=240 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p free-gpu -A rhayes1_lab_gpu --gres=gpu:1 --export=ALL --output=`pwd`/outrun --error=`pwd`/errrun ./RunSetup.sh

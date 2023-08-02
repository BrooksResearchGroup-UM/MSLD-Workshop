#! /bin/bash

# Set slurm options for cluster
# export SLURMOPTSMD="--time=1440 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p gpu --gres=gpu:1 --export=ALL"
# export SLURMOPTSPP="--time=1440 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p gpu --gres=gpu:1 --export=ALL"

# Run flattening
DEPEND=""
PID=`sbatch --parsable $SLURMOPTSMD $DEPEND ./runflat.sh`

# Run short production
DEPEND="--dependency=afterok:$PID"
NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do
  export step=61
  export p
  export nitt=5
  PID=`sbatch --parsable $SLURMOPTSMD --array=1-2%1 $DEPEND ./runprod.sh`
  NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
  COMMA=","
done

# Run postprocessing on short production to optimize biases
DEPEND=$NEWDEPEND
export i=61
export eqS=2
export S=10
export N=5
export skipE=1
PID=`sbatch --parsable $SLURMOPTSPP $DEPEND ./postprocess.sh`

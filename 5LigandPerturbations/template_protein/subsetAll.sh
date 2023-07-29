#! /bin/bash

# Set slurm options for cluster
export SLURMOPTSMD="--time=2-00:00:00 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p work --gres=gpu:1 --export=ALL"
export SLURMOPTSPP="--time=2-00:00:00 --ntasks=1 --tasks-per-node=1 --cpus-per-task=1 -p work --gres=gpu:1 --export=ALL"

# Run flattening
DEPEND=""
PID=`sbatch --parsable $SLURMOPTSMD $DEPEND ./runflat.sh`; echo $PID


# Run short production
DEPEND="--dependency=afterok:$PID"
NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do
  export step=114
  export p
  export nitt=5
  PID=`sbatch --parsable $SLURMOPTSMD --array=1-1%1 $DEPEND ./runprod.sh`; echo $PID
  NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
  COMMA=","
done

# Run postprocessing on short production to optimize biases
DEPEND=$NEWDEPEND
export i=114
export eqS=1
export S=5
export N=5
export skipE=1
PID=`sbatch --parsable $SLURMOPTSPP $DEPEND ./postprocess.sh`; echo $PID

# Run medium production
DEPEND="--dependency=afterok:$PID"
NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do
  export step=115
  export p
  export nitt=5
  PID=`sbatch --parsable $SLURMOPTSMD --array=1-5%1 $DEPEND ./runprod.sh`; echo $PID
  NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
  COMMA=","
done

# Run postprocessing on medium production to optimize biases
DEPEND=$NEWDEPEND
export i=115
export eqS=5
export S=25
export N=5
export skipE=1
PID=`sbatch --parsable $SLURMOPTSPP $DEPEND ./postprocess.sh`; echo $PID

# Run long production
DEPEND="--dependency=afterok:$PID"
NEWDEPEND="--dependency="
COMMA=""
for p in a b c d e
do
  export step=116
  export p
  export nitt=5
  PID=`sbatch --parsable $SLURMOPTSMD --array=1-10%1 $DEPEND ./runprod.sh`; echo $PID
  NEWDEPEND="${NEWDEPEND}${COMMA}afterok:$PID"
  COMMA=","
done

# Run postprocessing on long production to optimize biases and obtain results
DEPEND=$NEWDEPEND
export i=116
export eqS=10
export S=50
export N=5
export skipE=1
PID=`sbatch --parsable $SLURMOPTSPP $DEPEND ./postprocess.sh`; echo $PID

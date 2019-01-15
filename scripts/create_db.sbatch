#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --qos=normal
#SBATCH --output="create_db_%j.out"
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date # timestamp

# erinnishgrp@colostate.edu projects setup
PROJ_DIR=/projects/dcking@colostate.edu
source $PROJ_DIR/paths.bashrc

python3 create_db.py

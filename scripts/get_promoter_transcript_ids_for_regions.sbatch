#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:10:00
#SBATCH --qos=normal
NTHREADS=${SLURM_NTASKS} # passes --ntasks set above
echo "[$0] $SLURM_JOB_NAME $@" # log the command line
export TMPDIR=$SLURM_SCRATCH
export TMP=$TMPDIR
date # timestamp

# erinnishgrp@colostate.edu projects setup
PROJ_DIR=/projects/dcking@colostate.edu
source $PROJ_DIR/paths.bashrc

infile=$1
outfile=${infile/.bedlike/.gene_ids}
cmd="python3 $PRJ/WS263-repo/scripts/get_promoter_transcript_ids_for_regions.py $infile $PRJ/WS263/WS263.sqlite > $outfile"
echo $cmd
time eval $cmd

gene_id_file=${outfile}.filtered

#cmd2="grep -v NA$ $outfile | cut -f 12 | collapse , > $gene_id_file"
cmd2="grep -v NA$ $outfile | cut -f 12 | sed 's/Transcript://' > $gene_id_file"
echo "$cmd2"
eval $cmd2

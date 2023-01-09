#!/usr/bin/env bash

if [ -z "$path" ]
then
      CURR_DIR=$1
else
      CURR_DIR=$path
fi											# obtain current script directory

CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)					# obtain configuration file directory

source $CONFIG/general/config.sh							# load variable configuration file

ALIGNMENTS="$BASEDIR/analysis/alignments"
WD="$BASEDIR/analysis/gRNA_mods"
TEMP_DIR="$TEMPORARY/gRNA_mods"			


# parameters to edit
THREADS=12										# number of threads used for the scripts
PARALLEL_JOBS=5										# number of parallel jobs to run f5c
BASECALLING="guppy_initial"								# indicate basecalling version (since the initial basecalling for the different datasets has been performed with different versions of guppy, we just call it "guppy_initial")
condition="WT"										# condition of the samples to be processed


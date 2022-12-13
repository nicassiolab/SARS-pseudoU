#!/usr/bin/env bash

CURR_DIR=$(dirname "$(realpath "$0")")							# obtain current script directory
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)					# obtain configuration file directory

source $CONFIG/general/config.sh							# load variable configuration file

ALIGNMENTS="$BASEDIR/analysis/alignments"
WD="$BASEDIR/analysis/gRNA_mods"
TEMP_DIR="$TEMPORARY/gRNA_mods"			


# parameters to edit
THREADS=5
BASECALLING="guppy_initial"  #indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
condition="WT"


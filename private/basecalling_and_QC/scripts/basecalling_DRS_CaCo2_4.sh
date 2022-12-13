#!/bin/bash

# load variables
CURR_DIR=$(dirname "$(realpath "$0")")                                                  # obtain current script directory
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)                                    # obtain configuration file directory

source $CONFIG/general/config.sh                                                        # load variable configuration file
source $CURR_DIR/images.sh								# load images


WD="$TEMPORARY/basecalling/cov"
HAC_FILE="$FILES/rna_r9.4.1_70bps_hac.cfg"

# parameters to edit
THREADS=10
MIN_Q_SCORE=7
KIT="SQK-RNA002"
FLOWCELL="FLO-MIN106"

# basecalling with last guppy version of sample DRS_CaCo2_4, WT and KD 

mkdir -p $WD/c37_WT/fast5
cp $WT_raw/fast5_pass/*.fast5 $WD/c37_WT/fast5
cp $WT_raw/fast5_fail/*.fast5 $WD/c37_WT/fast5
$SINGC guppy_basecaller -i $WD/c37_WT/fast5 -s $WD --num_callers $THREADS -c $HAC_FILE --recursive --min_qscore $MIN_Q_SCORE --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna

mkdir -p $WD/c37_PUS7KD/fast5
cp $PUS7_KD_raw/fast5_pass/*.fast5 $WD/c37_PUS7KD/fast5
cp $PUS7_KD_raw/fast5_fail/*.fast5 $WD/c37_PUS7KD/fast5

$SINGC guppy_basecaller -i $WD/c37_PUS7KD/fast5 -s $WD --num_callers $THREADS -c $HAC_FILE --recursive --min_qscore $MIN_Q_SCORE --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna


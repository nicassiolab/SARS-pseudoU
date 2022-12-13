#!/bin/bash

# load variables
CURR_DIR=$(dirname "$(realpath "$0")")                                                  # obtain current script directory
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)                                    # obtain configuration file directory

source $CONFIG/general/config.sh                                                        # load variable configuration file
source $CURR_DIR/images.sh                                                             # load images


WD="$TEMPORARY/basecalling/cov"
HAC_FILE="$FILES/rna_r9.4.1_70bps_hac.cfg"

# parameters
THREADS=12
MIN_Q_SCORE=7
KIT="SQK-RNA002"
FLOWCELL="FLO-MIN106"

# basecalling with last guppy version of sample IVT 

mkdir -p $WD/IVT/fast5
cp $IVT-raw/*.fast5 $WD/IVT/fast5
$SINGC guppy_basecaller -i $WD/IVT/fast5 -s $WD --cpu_threads_per_caller 4 --num_callers $THREADS -c $HAC_FILE --recursive --min_qscore $MIN_Q_SCORE --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna



#!/bin/bash
set -e -o pipefail

# load variables
CURR_DIR=$(dirname "$(realpath "$0")")                                                  # obtain current script directory
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)                                    # obtain configuration file directory

source $CONFIG/general/config.sh                                                        # load variable configuration file
source $CURR_DIR/images.sh

THREADS=10
BASECALLING="guppy_initial"
WD="$BASEDIR/analysis/basecalling/pycoQC/$BASECALLING"



for condition in WT PUS7KD;do

	# take sample file for each condition and basecalling version
        SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"
	

	# read sample file 
        while IFS=$'\t' read sample fasta fast5 cell_line source; do
		SAMPLE="$sample"
                SAMPLE_FA="$fasta"
                SAMPLE_FAST5="$fast5"
                SAMPLE_CELL_LINE="$cell_line"
		mkdir -p $WD/$cell_line/"$condition"/
		
		# obtain guppy directory path and sequencing summary file
		GUPPY_DIR=${SAMPLE_FAST5%/*} 
	        SEQ_SUMMARY=$(find $GUPPY_DIR -name 'sequencing_summary*')
		
		# pycoQC quality controls
		if [ -f "$SEQ_SUMMARY" ]; then
			$PYCOQC pycoQC -f $SEQ_SUMMARY -o $WD/$cell_line/"$condition"/"$sample".html
		fi

	done < <(tail -n +2 $SAMPLE_FILE)

done



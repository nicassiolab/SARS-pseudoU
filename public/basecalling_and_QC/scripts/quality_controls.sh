#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis/basecalling/pycoQC"
FASTA="$BASEDIR/analysis/fasta"
IMG="$BASEDIR/img"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"
MOUNT_DIR="/hpcnfs"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
THREADS=10

# pull image
if [ ! -f "$IMG/pycoqc_2.5.2.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/pycoqc:2.5.2
fi
# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/pycoqc_2.5.2.sif"


# indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
BASECALLING="guppy_initial"
WD="$WD/$BASECALLING"



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
			$SINGC  pycoQC -f $SEQ_SUMMARY -o $WD/$cell_line/"$condition"/"$sample".html
		fi

	done < <(tail -n +2 $SAMPLE_FILE)

done



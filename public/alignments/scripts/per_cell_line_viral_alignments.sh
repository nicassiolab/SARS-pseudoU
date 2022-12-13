#!/bin/bash
set -e -o pipefail

# load variables from general configuration file
CURR_DIR=$(dirname "$(realpath "$0")")                                                  # obtain current script directory
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)                                    # obtain configuration file directory
source $CONFIG/general/config.sh

# load local configuration file
source $CURR_DIR/config.sh
# load images
source $CURR_DIR/images.sh

WD="$BASEDIR/analysis/alignments/$BASECALLING/per_cell_line"
WD_MODS="$BASEDIR/analysis/sgRNAs_mods_detection/$BASECALLING/per_cell_line"
FASTA="$BASEDIR/analysis/fasta/$BASECALLING"
FAST5="$BASEDIR/analysis/fast5/per_cell_line"


# take sample file for basecalling version
SAMPLE_FILE="${FILES}/${condition_per_cell_line}_samples_${BASECALLING}.txt"

# read sample file
for selected_cell_line in CaCo2 CaLu3 VeroE6; do
	mkdir -p $FASTA/per_cell_line/$condition_per_cell_line
	find $FASTA/$selected_cell_line -type f -name '*.fa' -exec cat {} \; > $FASTA/per_cell_line/$condition_per_cell_line/"$selected_cell_line".fa

	# align fasta to the reference transcriptome
	mkdir -p $WD/$condition_per_cell_line/alignments_to_assembly/alignments_backup/
	$SINGC minimap2 -t $THREADS -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $FASTA/per_cell_line/$condition_per_cell_line/"$selected_cell_line".fa > $WD/"$condition_per_cell_line"/alignments_to_assembly/"$selected_cell_line".sam
	$SINGC samtools view -h $WD/"$condition_per_cell_line"/alignments_to_assembly/"$selected_cell_line".sam > $WD/"$condition_per_cell_line"/alignments_to_assembly/alignments_backup/"$selected_cell_line".bam
	$SINGC samtools view -h -F 2324 -Sb $WD/"$condition_per_cell_line"/alignments_to_assembly/"$selected_cell_line".sam | $SINGC samtools sort > $WD/"$condition_per_cell_line"/alignments_to_assembly/"$selected_cell_line"_nanocompore.bam
	$SINGC samtools index $WD/"$condition_per_cell_line"/alignments_to_assembly/"$selected_cell_line"_nanocompore.bam
	rm $WD/"$condition_per_cell_line"/alignments_to_assembly/"$selected_cell_line".sam

done



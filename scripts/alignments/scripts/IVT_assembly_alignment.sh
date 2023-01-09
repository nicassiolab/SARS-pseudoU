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

WD="$BASEDIR/analysis/alignments"
SAMPLE_FILE="${FILES}/WT_samples_${BASECALLING}.txt"

while IFS=$'\t' read sample fasta fast5 cell_line source; do

	mkdir -p $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/bam
	mkdir -p $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/sorted
	$SINGC minimap2 -x map-ont -t $THREADS -a $TRANSCRIPTOME_ASSEMBLY $fasta > $WD/$BASECALLING/IVT/WT/alignments_to_assembly/IVT.sam
	$SINGC samtools view -h -Sb $WD/$BASECALLING/IVT/WT/alignments_to_assembly/IVT.sam > $WD/$BASECALLING/IVT/WT/alignments_to_assembly/IVT.bam
	$SINGC samtools view -h -F 2068 $WD/$BASECALLING/IVT/WT/alignments_to_assembly/IVT.bam > $WD/$BASECALLING/IVT/WT/alignments_to_assembly/IVT_filtered.bam
	$SINGC samtools view -F 2068 $WD/$BASECALLING/IVT/WT/alignments_to_assembly/IVT.sam > $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/IVT_filtered.sam
	cd $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/
	awk '{print>$3".sam"}' $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/IVT_filtered.sam
	rm $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/IVT_filtered.sam
	for file in $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/*.sam ; do mv $file ${file//|/_} ; done
	for filename in $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/*.sam; do
		name=${filename##*/}
		base=${name%.sam}
		$SINGC samtools view -h -t $TRANSCRIPTOME_ASSEMBLY_FAI -Sb $filename > $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/bam/"$base".bam
		$SINGC samtools sort $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/bam/"$base".bam > $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/sorted/"$base"_sorted.bam
		$SINGC samtools index $WD/$BASECALLING/IVT/WT/alignments_to_assembly/splitted/sorted/"$base"_sorted.bam
	done

done < <(grep "IVT" $SAMPLE_FILE)


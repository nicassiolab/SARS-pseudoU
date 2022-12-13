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


# directories
WD="$BASEDIR/analysis/alignments/$BASECALLING"
FASTA="$BASEDIR/analysis/fasta/$BASECALLING"



for condition in $SAMPLE_CONDITION ;do

	# take sample file for each condition and basecalling version
        SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"
	

	# read sample file 
        while IFS=$'\t' read sample fasta fast5 cell_line source; do
		SAMPLE="$sample"
                SAMPLE_FA="$fasta"
                SAMPLE_FAST5="$fast5"
                SAMPLE_CELL_LINE="$cell_line"
		mkdir -p $FASTA/$cell_line/"$condition"/

		# align fasta to the reference human  genome
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_human_genome/alignments_backup/
	        $SINGC minimap2 $HG_GENOME_PARAM -t $THREADS $GENOME_FA $FASTA/$cell_line/"$condition"/"$sample".fa > $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample".sam
        	$SINGC samtools view -h $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample".sam > $WD/$cell_line/"$condition"/alignments_to_human_genome/alignments_backup/"$sample".bam
        	$SINGC samtools view -F 2324 -Sb $WD/$cell_line/"$condition"/alignments_to_human_genome/alignments_backup/"$sample".bam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample"_sorted.bam
        	rm $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample".sam
        	$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample"_sorted.bam
        	$SINGC bedtools bamtobed -bed12 -i $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample"_sorted.bam > $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample"_sorted.bed

		# align fasta to the reference transcriptome
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/alignments_backup/
		$SINGC minimap2 $TRANSCRIPTOME_PARAM -t $THREADS $HG_DATA/transcriptome_fasta.fa $FASTA/$cell_line/"$condition"/"$sample".fa > $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample".sam
		$SINGC samtools view -h $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample".sam > $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/alignments_backup/"$sample".bam
		$SINGC samtools view -F 2324 -Sb $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/alignments_backup/"$sample".bam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocompore.bam
		$SINGC samtools view -F 2068 -Sb $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample".sam |  $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocount.bam
		$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocount.bam
		$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocompore.bam

	done < <(tail -n +2 $SAMPLE_FILE)

done



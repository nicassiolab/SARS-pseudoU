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


	        # transform fastq file into fasta file and copy to directory
	        if [ -d $SAMPLE_FA ]; then
        	        cat $SAMPLE_FA/*.fastq > $FASTA/$cell_line/"$condition"/"$sample".fastq
                        awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' $FASTA/$cell_line/"$condition"/"$sample".fastq > $FASTA/$cell_line/"$condition"/"$sample".fa
                        rm $FASTA/$cell_line/"$condition"/"$sample".fastq
        	elif [ -f $SAMPLE_FA ]; then
                	cp $SAMPLE_FA $FASTA/$cell_line/"$condition"/"$sample".fa
        	fi

		# align fasta to the reference genome
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_genome/alignments_backup/
	        $SINGC minimap2 $GENOME_PARAM -t $THREADS $GENOME_FA $FASTA/$cell_line/"$condition"/"$sample".fa > $WD/$cell_line/"$condition"/alignments_to_genome/"$sample".sam
        	$SINGC samtools view -h $WD/$cell_line/"$condition"/alignments_to_genome/"$sample".sam > $WD/$cell_line/"$condition"/alignments_to_genome/alignments_backup/"$sample".bam
        	$SINGC samtools view -F 2324 -Sb $WD/$cell_line/"$condition"/alignments_to_genome/alignments_backup/"$sample".bam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bam
        	rm $WD/$cell_line/"$condition"/alignments_to_genome/"$sample".sam
        	$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bam
        	$SINGC bedtools bamtobed -bed12 -i $WD/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bam > $WD/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bed

		# align fasta to the reference transcriptome
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_assembly/alignments_backup/
		$SINGC minimap2 -t $THREADS $TRANSCRIPTOME_PARAM $TRANSCRIPTOME_ASSEMBLY $FASTA/$cell_line/"$condition"/"$sample".fa > $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam
		$SINGC samtools view -h $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam > $WD/$cell_line/"$condition"/alignments_to_assembly/alignments_backup/"$sample".bam
		$SINGC samtools view -h -F 2068 -Sb $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocount.bam
		$SINGC samtools view -h -F 2324 -Sb $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocompore.bam
		$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocount.bam
		$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocompore.bam
		rm $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam


	done < <(tail -n +2 $SAMPLE_FILE)

done



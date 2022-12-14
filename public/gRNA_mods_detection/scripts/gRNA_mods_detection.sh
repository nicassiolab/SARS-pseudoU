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


# assign working directories to a variable  
WD="$WD/$BASECALLING/nanocompore"
FASTA="$FASTA/$BASECALLING"

# take sample file for each condition and basecalling version
SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"


# select all genome-length reads from alignments
while IFS=$'\t' read sample fasta fast5 cell_line source; do
	cat $FASTA/$cell_line/"$condition"/"$sample".fa >> $TEMP_DIR/all_samples.fa
	mkdir -p $TEMP_DIR
        cp $sample_fast5/*.fast5 $TEMP_DIR/fast5       # this line will be used for eventalign

	mkdir -p $WD/$cell_line/"$sample"_gRNAs.txt	                
	awk '{if($2<=45 && $3>=29850 && $10==1) print $4}' $ALIGNMENTS/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bed > $WD/$cell_line/"$sample"_gRNAs.txt
	$SINGC python3 $SCRIPTDIR/extract_reads_bam.py $ALIGNMENTS/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bam $WD/$cell_line/"$sample"_gRNAs.bam $WD/$cell_line/"$sample"_gRNAs.txt
	rm $WD/$cell_line/"$sample"_gRNAs.txt

done < <(awk '$1!="IVT"' $SAMPLE_FILE | tail -n +2)


mkdir -p $WD/all_cell_lines/alignments_to_genome
bam_list=$(find $WD/*/ -name "*_gRNAs.bam")
$SINGC samtools merge $WD/all_cell_lines/alignments_to_genome/gRNAs.bam $bam_list
$SINGC samtools view -h -F 2324 -Sb $WD/all_cell_lines/alignments_to_genome/gRNAs.bam | $SINGC samtools sort > $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam
$SINGC samtools index $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam


# eventalign reads
mkdir -p $WD/all_cell_lines/eventalign/collapse
$F5C f5c index -t $THREADS --iop $PARALLEL_JOBS -d $TEMP_DIR/fast5 $TEMP_DIR/all_samples.fa
$F5C f5c eventalign --rna --min-mapq 0 -t $THREADS -r $TEMP_DIR/all_samples.fa -b $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam --g $GENOME_FA --samples --print-read-names --scale-events --iop 5  | $NANOPOLISHCOMP NanopolishComp Eventalign_collapse -t $THREADS -o $WD/all_cell_lines/eventalign/collapse


# IVT processing
mkdir -p $WD/IVT
$F5C f5c eventalign --rna --min-mapq 0 -t $THREADS -r $ANALYSIS/fasta/IVT/IVT.fa -b $ANALYSIS/alignments/IVT/alignments_to_genome/IVT_sorted.bam --g $GENOME_FA --samples --print-read-names --scale-events --iop $PARALLEL_JOBS  | $NANOPOLISHCOMP NanopolishComp Eventalign_collapse -t $THREADS -o $WD/IVT/eventalign/collapse


# nanocompore
$NANOCOMPORE nanocompore sampcomp --file_list1 $WD/all_cell_lines/eventalign/collapse/out_eventalign_collapse.tsv --file_list2 $WD/IVT/eventalign/collapse/out_eventalign_collapse.tsv --label1 WT --label2 IVT --fasta $GENOME_FA --outpath $WD/all_cell_lines/nanocompore --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads $THREADS



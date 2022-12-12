#!/bin/bash
set -e -o pipefail


##################################################################################### PUS7KD vs all WT samples ############################################################################################

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"
ALIGNMENTS="$BASEDIR/analysis/alignments"
WD="$BASEDIR/analysis/gRNA_mods"


BASECALLING="guppy_v601"
condition="PUS7KD"
WD="$WD/$BASECALLING/nanocompore/KD_vs_all"
FASTA="$FASTA/$BASECALLING"
SAMPLE="DRS_CaCo2_4"

# take sample file for each condition and basecalling version
SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"

# select all genome-length reads from alignments
while IFS=$'\t' read sample fasta fast5 cell_line source; do

        mkdir -p $WD/$condition/
        awk '{if($2<=45 && $3>=29850 && $10==1) print $4}' $ALIGNMENTS/$BASECALLING/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bed > $WD/$condition/"$sample"_gRNAs.txt
        awk '{split($11,a,",");sum = 0;for( i = 1; i <= length( a ); i++ ) sum += a[i]; print sum}' $ALIGNMENTS/$BASECALLING/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bed | awk '$1>=28000' | wc -l > $WD/$condition/"$sample"_gRNAs_number.txt            # evaluate number of reads with length > 28000 nt (full-length genomic reads)

done < <(grep $SAMPLE $SAMPLE_FILE)


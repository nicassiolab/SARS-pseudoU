#!/bin/bash

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
WD="${BASEDIR}/analysis/per_cell_line/shared"
DATADIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/shared/tracks_for_overlap"

mkdir -p $WD/in_2
mkdir -p $WD/in_3

###	overlap over 3 databases

source activate $HOME/miniconda/envs/nanopore_pipeline_env/
for dir in $DATADIR/*; do
	echo $dir
	bedtools merge -s -i $dir/calu3.bed > $dir/calu3_merged.bed
	bedtools merge -s -i $dir/vero.bed > $dir/vero_merged.bed
	bedtools merge -s -i $dir/caco2.bed > $dir/caco2_merged.bed
	intersectBed -a $dir/calu3_merged.bed -b $dir/caco2_merged.bed > $dir/cc.bed
	intersectBed -a $dir/cc.bed -b $dir/vero_merged.bed > $dir/common.bed
	bedtools merge -i $dir/common.bed > $dir/common_merged.bed3
	#rm $dir/*.bed

done



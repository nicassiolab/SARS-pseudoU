#!/bin/bash
set -e -o pipefail


WD="$BASEDIR/analysis/sgRNAs_mods_detection/$BASECALLING"
FASTA="$BASEDIR/analysis/fasta/$BASECALLING"
ALIGNMENTS="$BASEDIR/analysis/alignments/$BASECALLING"


for condition in WT PUS7KD;do

	# take sample file for each condition and basecalling version
        SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"
	
	# read sample file 
        while IFS=$'\t' read sample fasta fast5 cell_line source; do
		SAMPLE="$sample"
                SAMPLE_FA="$fasta"
                SAMPLE_FAST5="$fast5"
                SAMPLE_CELL_LINE="$cell_line"
		mkdir -p $FASTA/$cell_line/"$condition"/


		$NANOCOMPORE f5c index -t $THREADS --iop 5 -d $SAMPLE_FAST5 $FASTA/$cell_line/"$condition"/"$sample".fa
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/eventalign/
		$NANOCOMPORE sh -c "f5c eventalign --rna --min-mapq 0 -t $THREADS -r $FASTA/$cell_line/"$condition"/"$sample".fa -b $ALIGNMENTS/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocompore.bam --g $HG_DATA/transcriptome_fasta.fa --samples --print-read-names --disable-cuda=yes --scale-events --iop $PARALLEL_JOBS | nanocompore eventalign_collapse -o $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/eventalign/collapse"



	done < <(tail -n +2 $SAMPLE_FILE)

done





# indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
BASECALLING="guppy_v601"
WD="$WD/$BASECALLING"
condition_1="WT"
condition_2="PUS7KD"
cell_line="CaCo2"

$NANOCOMPORE nanocompore sampcomp --file_list1 $WD/$cell_line/$condition_1/alignments_to_human_transcriptome/eventalign/collapse/out_eventalign_collapse.tsv --file_list2 $WD/$cell_line/$condition_2/alignments_to_human_transcriptome/eventalign/collapse/out_eventalign_collapse.tsv  --label1 $condition_1 --label2 $condition_2 --fasta $NANOCOMP_FA --outpath $WD/$cell_line/nanocompore/alignments_to_human_transcriptome/"$condition_1"_vs_"$condition_2"/ --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads $THREADS --bed $NANOCOMP_BED



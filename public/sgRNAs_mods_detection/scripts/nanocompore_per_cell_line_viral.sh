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

WD="$BASEDIR/analysis/sgRNAs_mods_detection/$BASECALLING"
FASTA="$BASEDIR/analysis/fasta/$BASECALLING"


# take sample file for basecalling version
SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"


# read sample file
for selected_cell_line in CaCo2 VeroE6 CaLu3; do
	
	# copy all raw fast5 in a directory for each cell line
	mkdir -p $WD/fast5_temporary/$selected_cell_line
		
	while IFS=$'\t' read sample fasta fast5 cell_line source; do
		cp $fast5/*.fast5 $WD/fast5_temporary/$selected_cell_line
	done < <(grep $selected_cell_line $SAMPLE_FILE)

	# index and eventalign files for each cell line
	$NANOCOMPORE f5c index -d $WD/fast5_temporary/$selected_cell_line $FASTA/$condition/"$selected_cell_line".fa
	mkdir -p $WD/$selected_cell_line/$condition/eventalign
	$NANOCOMPORE sh -c "f5c eventalign --rna --min-mapq 0 -t $THREADS -r $FASTA/$condition/"$selected_cell_line".fa -b $WD/alignments/"$condition"/alignments_to_assembly/"$selected_cell_line"_nanocompore.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events --disable-cuda=yes --iop $PARALLEL_JOBS |  nanocompore eventalign_collapse -o $WD/$selected_cell_line/$condition/eventalign/collapse"

	# compare eventalign files for each cell line to eventalign files for IVT, sgRNA per sgRNA
	for filename in $WD/IVT/eventalign*; do
		base=${filename##*/eventalign_}
		$NANOCOMPORE nanocompore sampcomp --file_list1 $WD/$selected_cell_line/$condition/eventalign/collapse/out_eventalign_collapse.tsv  --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 $condition --label2 IVT --fasta $TRANSCRIPTOME_ASSEMBLY --outpath  $WD/$selected_cell_line/nanocompore/sampcomp/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads $THREADS --bed $NANOCOMP_BED

	done 

done


# generate violin plots 
$SINGC python3 $CURR_DIR/pseudoU_signal_plots.py --transcriptome_assembly $TRANSCRIPTOME_ASSEMBLY --basecalling_version $BASECALLING --base_directory $BASEDIR



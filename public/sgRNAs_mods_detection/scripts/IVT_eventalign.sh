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
SAMPLE_FILE="${FILES}/WT_samples_${BASECALLING}.txt"

# index and eventalign files for IVT

while IFS=$'\t' read sample fasta fast5 cell_line source; do

	$NANOCOMPORE f5c index -d $fast5 $fasta
	for filename in $BASEDIR/analysis/alignments/IVT/WT/alignments_to_assembly/splitted/sorted/*.bam; do
		name=${filename##*/}
                base=${name%.bam}
	
		mkdir -p $WD/IVT/eventalign_"$base"
		$NANOCOMPORE sh -c "f5c eventalign --rna --min-mapq 0 -t $THREADS -r $fasta -b $BASEDIR/analysis/alignments/IVT/WT/alignments_to_assembly/splitted/sorted/"$base".bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events --disable-cuda=yes --iop $PARALLEL_JOBS |  nanocompore eventalign_collapse -o $WD/IVT/eventalign_"$base""
	done

done < <(grep "IVT" $SAMPLE_FILE)



#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis/sgRNAs_mods_detection"
ALIGNMENTS="$BASEDIR/analysis/alignments"
FASTA="$BASEDIR/analysis/fasta"
IMG="$BASEDIR/img"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"
MOUNT_DIR="/hpcnfs/scratch"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
THREADS=10


# pull image for NanopolishComp
if [ ! -f "$IMG/nanocompore_latest.sif" ]; then
        cd $IMG
        singularity pull docker://adrienleger/nanocompore:latest
fi
# singularity command
NANOPOLISHCOMP="singularity exec -B $MOUNT_DIR $IMG/nanocompore_latest.sif"

# pull image for nanocompore
if [ ! -f "$IMG/nanocompore_v1.0.4.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/nanocompore:v1.0.4
fi
# singularity command
NANOCOMPORE="singularity exec -B $MOUNT_DIR $IMG/nanocompore_v1.0.4.sif"

# pull image for f5c
if [ ! -f "$IMG/f5c_v0.6.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/f5c:v0.6
fi
# singularity command
F5C="singularity exec -B $MOUNT_DIR $IMG/f5c_v0.6.sif"


# indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
BASECALLING="guppy_v601"
WD="$WD/$BASECALLING"
FASTA="$FASTA/$BASECALLING"


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


		$F5C f5c index -d $SAMPLE_FAST5 $FASTA/$cell_line/"$condition"/"$sample".fa
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_assembly/eventalign/
		$F5C f5c eventalign --rna --min-mapq 0 -t $THREADS -r $SAMPLE_FAST5 $FASTA/$cell_line/"$condition"/"$sample".fa -b $ALIGNMENTS/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocompore.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events | $NANOPOLISHCOMP NanopolishComp Eventalign_collapse -o $WD/$cell_line/"$condition"/alignments_to_assembly/eventalign/



	done < <(tail -n +2 $SAMPLE_FILE)

done



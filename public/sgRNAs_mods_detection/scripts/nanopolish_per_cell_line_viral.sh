#!/bin/bash
set -e -o pipefail

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis/sgRNAs_mods_detection"
FASTA="$BASEDIR/analysis/fasta"
IMG="$BASEDIR/img"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"
MOUNT_DIR="/hpcnfs/scratch"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/consensus_extracted.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"
THREADS=10

# pull image for nanocompore
if [ ! -f "$IMG/nanocompore_v1.0.4.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/nanocompore:v1.0.4
fi
# singularity command
NANOCOMPORE="singularity exec -B $MOUNT_DIR $IMG/nanocompore_v1.0.4.sif"

# indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
BASECALLING="guppy_v601"
WD="$WD/$BASECALLING/per_cell_line"
FASTA="$FASTA/$BASECALLING/per_cell_line"
condition="WT"

# take sample file for basecalling version
SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"


# read sample file
for selected_cell_line in CaCo2 VeroE6 CaLu3; do
	
	mkdir -p $WD/fast5_temporary/$selected_cell_line
	
	while IFS=$'\t' read sample fasta fast5 cell_line source; do
		cp $fast5/*.fast5 $WD/fast5_temporary/$selected_cell_line
	done < <(grep $selected_cell_line $SAMPLE_FILE)

	$NANOCOMPORE f5c index -d $WD/fast5_temporary/$selected_cell_line $FASTA/$condition/"$selected_cell_line".fa
	mkdir -p $WD/$selected_cell_line/$condition/eventalign
	$NANOCOMPORE sh -c "f5c eventalign --rna --min-mapq 0 -t $THREADS -r $FASTA/$condition/"$selected_cell_line".fa -b $WD/alignments/"$condition"/alignments_to_assembly/"$selected_cell_line"_nanocompore.bam --g $TRANSCRIPTOME_ASSEMBLY --samples --print-read-names --scale-events --disable-cuda=yes --iop 5 |  nanocompore eventalign_collapse -o $WD/$selected_cell_line/$condition/eventalign/collapse"

done





#!/bin/bash
set -e -o pipefail

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis/alignments"
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


# pull image
if [ ! -f "$IMG/nrceq_pipeline_latest.sif" ]; then
        cd $IMG
        singularity pull docker://cugolini/nrceq_pipeline:latest
fi


# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/nrceq_pipeline_latest.sif"

# indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
BASECALLING="guppy_v601"
WD="$WD/$BASECALLING/per_cell_line"
FASTA="$FASTA/$BASECALLING"
condition="WT"


# take sample file for basecalling version
SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"

# read sample file
for selected_cell_line in CaCo2 IVT VeroE6; do
	mkdir -p $FASTA/per_cell_line/$condition
	find $FASTA/$selected_cell_line -type f -name '*.fa' -exec cat {} \; > $FASTA/per_cell_line/$condition/"$selected_cell_line".fa

	# align fasta to the reference transcriptome
	mkdir -p $WD/$condition/alignments_to_assembly/alignments_backup/
	$SINGC minimap2 -t $THREADS -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $FASTA/per_cell_line/$condition/"$selected_cell_line".fa > $WD/"$condition"/alignments_to_assembly/"$selected_cell_line".sam
	$SINGC samtools view -h $WD/"$condition"/alignments_to_assembly/"$selected_cell_line".sam > $WD/"$condition"/alignments_to_assembly/alignments_backup/"$selected_cell_line".bam
	$SINGC samtools view -h -F 2324 -Sb $WD/"$condition"/alignments_to_assembly/"$selected_cell_line".sam | $SINGC samtools sort > $WD/"$condition"/alignments_to_assembly/"$selected_cell_line"_nanocompore.bam
	$SINGC samtools index $WD/"$condition"/alignments_to_assembly/"$selected_cell_line"_nanocompore.bam
	rm $WD/"$condition"/alignments_to_assembly/"$selected_cell_line".sam

done




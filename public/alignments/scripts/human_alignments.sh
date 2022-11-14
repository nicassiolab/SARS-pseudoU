#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis/alignments"
FASTA="$BASEDIR/analysis/fasta"
IMG="$BASEDIR/img"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"
MOUNT_DIR="/hpcnfs/scratch"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/Datasets/HUMAN_REFERENCE/Homo_sapiens.GRCh38.dna_sm.toplevel.fa"
GENOME_PARAM="-ax splice -uf -k14 -I100g"
TRANSCRIPTOME_FA="/hpcnfs/scratch/TSSM/cugolini/Datasets/HUMAN_REFERENCE/transcriptome_fasta.fa"
TRANSCRIPTOME_PARAM="-ax map-ont -p 0 -N 10" 
THREADS=10


# pull image
if [ ! -f "$IMG/nrceq_pipeline_latest.sif" ]; then
        cd $IMG
        singularity pull docker://cugolini/nrceq_pipeline:latest
fi
# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/nrceq_pipeline_latest.sif"

# pull image of samtools 1.9 for f5c
if [ ! -f "$IMG/samtools_1.9--76b9270.sif" ]; then
        cd $IMG
        singularity pull docker://nanozoo/samtools:1.9--76b9270
fi
SINGSAM="singularity exec -B $MOUNT_DIR $IMG/samtools_1.9--76b9270.sif"

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

		# align fasta to the reference human  genome
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_human_genome/alignments_backup/
	        $SINGC minimap2 $GENOME_PARAM -t $THREADS $GENOME_FA $FASTA/$cell_line/"$condition"/"$sample".fa > $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample".sam
        	$SINGC samtools view -h $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample".sam > $WD/$cell_line/"$condition"/alignments_to_human_genome/alignments_backup/"$sample".bam
        	$SINGC samtools view -F 2324 -Sb $WD/$cell_line/"$condition"/alignments_to_human_genome/alignments_backup/"$sample".bam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample"_sorted.bam
        	rm $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample".sam
        	$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample"_sorted.bam
        	$SINGC bedtools bamtobed -bed12 -i $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample"_sorted.bam > $WD/$cell_line/"$condition"/alignments_to_human_genome/"$sample"_sorted.bed

		# align fasta to the reference transcriptome
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/alignments_backup/
		$SINGC minimap2 $TRANSCRIPTOME_PARAM -t $THREADS $TRANSCRIPTOME_FA $FASTA/$cell_line/"$condition"/"$sample".fa > $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample".sam
		$SINGC samtools view -h $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample".sam > $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/alignments_backup/"$sample".bam
		$SINGC samtools view -F 2324 -Sb $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/alignments_backup/"$sample".bam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocompore.bam
		$SINGC samtools view -F 2068 -Sb $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample".sam |  $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocount.bam
		$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocount.bam
		$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_human_transcriptome/"$sample"_nanocompore.bam

	done < <(tail -n +2 $SAMPLE_FILE)

done



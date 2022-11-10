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
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
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

	        # transform fastq file into fasta file and copy to directory
	        if [ -d $SAMPLE_FA ]; then
        	        cat $SAMPLE_FA/*.fastq > $FASTA/$cell_line/"$condition"/"$sample".fastq
                        awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' $FASTA/$cell_line/"$condition"/"$sample".fastq > $FASTA/$cell_line/"$condition"/"$sample".fa
                        rm $FASTA/$cell_line/"$condition"/"$sample".fastq
        	elif [ -f $SAMPLE_FA ]; then
                	cp $SAMPLE_FA $FASTA/$cell_line/"$condition"/"$sample".fa
        	fi

		# align fasta to the reference genome
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_genome/alignments_backup/
	        $SINGC minimap2 $GENOME_PARAM -t $THREADS $GENOME_FA $FASTA/$cell_line/"$condition"/"$sample".fa > $WD/$cell_line/"$condition"/alignments_to_genome/"$sample".sam
        	$SINGC samtools view -h $WD/$cell_line/"$condition"/alignments_to_genome/"$sample".sam > $WD/$cell_line/"$condition"/alignments_to_genome/alignments_backup/"$sample".bam
        	$SINGC samtools view -F 2324 -Sb $WD/$cell_line/"$condition"/alignments_to_genome/alignments_backup/"$sample".bam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bam
        	rm $WD/$cell_line/"$condition"/alignments_to_genome/"$sample".sam
        	$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bam
        	$SINGC bedtools bamtobed -bed12 -i $WD/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bam > $WD/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bed

		# align fasta to the reference transcriptome
		mkdir -p $WD/$cell_line/"$condition"/alignments_to_assembly/alignments_backup/
		$SINGC minimap2 -t $THREADS -ax map-ont -p 0 -N 10 $TRANSCRIPTOME_ASSEMBLY $FASTA/$cell_line/"$condition"/"$sample".fa > $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam
		$SINGC samtools view -h $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam > $WD/$cell_line/"$condition"/alignments_to_assembly/alignments_backup/"$sample".bam
		$SINGC samtools view -h -F 2068 -Sb $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocount.bam
		$SINGC samtools view -h -F 2324 -Sb $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam | $SINGC samtools sort > $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocompore.bam
		$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocount.bam
		$SINGC samtools index $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample"_nanocompore.bam
		rm $WD/$cell_line/"$condition"/alignments_to_assembly/"$sample".sam

	done < <(tail -n +2 $SAMPLE_FILE)

done



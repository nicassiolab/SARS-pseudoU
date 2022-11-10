#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
ANALYSIS="$BASEDIR/analysis"
SCRIPTDIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/public/gRNA_mods_detection/scripts"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"
WD="$BASEDIR/analysis/gRNA_mods"
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
THREADS=10
TEMP_DIR="/hpcnfs/scratch/temporary/camilla_FN/gRNA_mods"


# pull image
if [ ! -f "$IMG/nrceq_pipeline_latest.sif" ]; then
        cd $IMG
        singularity pull docker://cugolini/nrceq_pipeline:latest
fi

# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/nrceq_pipeline_latest.sif"



# indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
BASECALLING="guppy_initial"
condition="WT"
WD="$WD/$BASECALLING"
FASTA="$FASTA/$BASECALLING"

# take sample file for each condition and basecalling version
SAMPLE_FILE="${FILES}/${condition}_samples_${BASECALLING}.txt"


# select all genome-length reads from alignments
while IFS=$'\t' read sample fasta fast5 cell_line source; do
	cat $FASTA/$cell_line/"$condition"/"$sample".fa >> $TEMP_DIR/all_samples.fa
        #cp $sample_fast5/*.fast5 $TEMP_DIR/fast5       # this line will be used for eventalign

	mkdir -p $WD/$cell_line/"$sample"_gRNAs.txt	                
	awk '{if($2<=45 && $3>=29850 && $10==1) print $4}' $ALIGNMENTS/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bed > $WD/$cell_line/"$sample"_gRNAs.txt
	$SINGC python3 $SCRIPTDIR/extract_reads_bam.py $ALIGNMENTS/$cell_line/"$condition"/alignments_to_genome/"$sample"_sorted.bam $WD/$cell_line/"$sample"_gRNAs.bam $WD/$cell_line/"$sample"_gRNAs.txt
	rm $WD/$cell_line/"$sample"_gRNAs.txt

done < <(awk '$1!="IVT"' $SAMPLE_FILE | tail -n +2)


mkdir -p $WD/all_cell_lines/alignments_to_genome
bam_list=$(find $WD/*/ -name "*_gRNAs.bam")
$SINGC samtools merge $WD/all_cell_lines/alignments_to_genome/gRNAs.bam $bam_list
$SINGC samtools view -h -F 2324 -Sb $WD/all_cell_lines/alignments_to_genome/gRNAs.bam | $SINGC samtools sort > $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam
$SINGC samtools index $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam


# eventalign reads
mkdir -p $WD/all_cell_lines/eventalign/collapse
$F5C f5c index -t $THREADS --iop 5 -d $TEMP_DIR/fast5 $TEMP_DIR/all_samples.fa
$F5C f5c eventalign --rna --min-mapq 0 -t $THREADS -r $TEMP_DIR/all_samples.fa -b $WD/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam --g $GENOME_FA --samples --print-read-names --scale-events --iop 5  | $NANOPOLISHCOMP NanopolishComp Eventalign_collapse -t $THREADS -o $WD/all_cell_lines/eventalign/collapse


# IVT processing
mkdir -p $WD/IVT
$F5C f5c eventalign --rna --min-mapq 0 -t $THREADS -r $ANALYSIS/fasta/IVT/IVT.fa -b $ANALYSIS/alignments/IVT/alignments_to_genome/IVT_sorted.bam --g $GENOME_FA --samples --print-r
ead-names --scale-events --iop 5  | $NANOPOLISHCOMP NanopolishComp Eventalign_collapse -t $THREADS -o $WD/IVT/eventalign/collapse


# nanocompore
$NANOCOMPORE nanocompore sampcomp --file_list1 $WD/all_cell_lines/eventalign/collapse/out_eventalign_collapse.tsv --file_list2 $WD/IVT/eventalign/collapse/out_eventalign_collapse.tsv --label1 WT --label2 IVT --fasta $GENOME_FA --outpath $WD/all_cell_lines/nanocompore --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads $THREADS

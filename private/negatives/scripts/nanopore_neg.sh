#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov/negatives"
DATA="$BASEDIR/data"
FILES="$BASEDIR/files"
WD="$BASEDIR/analysis"
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch/"
THREADS=8
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
NEGATIVES="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210423_1316_X3_FAL83184_52a41d7d/S33624_Calu3-CoV2-infected"
OLIGO_ONT_ADAPT="$FILES/ont_adapter.fa"


# pull image
if [ ! -f "$IMG/nrceq_pipeline_latest.sif" ]; then
        cd $IMG
        singularity pull docker://cugolini/nrceq_pipeline:latest
fi

# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/nrceq_pipeline_latest.sif"


#### analysis of negative-sense transcripts obtained by specific adapter at the genomic 5'

# pycoQC
if [ ! -f "$IMG/pycoqc_2.5.2.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/pycoqc
fi


mkdir -p $WD/negatives/pycoqc
singularity exec -B $MOUNT_DIR $IMG/pycoqc_2.5.2.sif pycoQC -f $NEGATIVES/sequencing_summary_FAL83184_c1102ce7.txt -o $WD/negatives/pycoqc/sequencing_summary_NEGATIVES.html

# read analysis (produce tracks)
mkdir -p $WD/negatives/fastq
mkdir -p $WD/negatives/alignments_to_genome/alignments_backup
mkdir -p $WD/negatives/alignments_to_genome/unmapped
cat $NEGATIVES/fastq_pass/*.fastq > $WD/negatives/fastq/pass.fastq
$SINGC minimap2 -t $THREADS $GENOME_FA $WD/negatives/fastq/pass.fastq > $WD/negatives/alignments_to_genome/negatives.sam
$SINGC samtools view -Sb $WD/negatives/alignments_to_genome/negatives.sam > $WD/negatives/alignments_to_genome/alignments_backup/negatives.bam
$SINGC samtools view -h -F 4 $WD/negatives/alignments_to_genome/alignments_backup/negatives.bam | $SINGC samtools sort > $WD/negatives/alignments_to_genome/negatives_filt_sort.bam
$SINGC samtools index $WD/negatives/alignments_to_genome/negatives_filt_sort.bam
rm $WD/negatives/alignments_to_genome/negatives.sam $WD/negatives/fastq/pass.fastq
$SINGC samtools view -h -f 16 $WD/negatives/alignments_to_genome/alignments_backup/negatives.bam | $SINGC samtools sort > $WD/negatives/alignments_to_genome/negatives_minus_filt_sort.bam
$SINGC samtools index $WD/negatives/alignments_to_genome/negatives_minus_filt_sort.bam

# analyse unmapped reads
UNMAPPED_DIR="$WD/negatives/alignments_to_genome/unmapped"

# parasail
if [ ! -f "$IMG/parasail_2.4.3.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/parasail:2.4.3 
fi


$SINGC samtools view -h -f 4 $WD/negatives/alignments_to_genome/alignments_backup/negatives.bam | $SINGC samtools sort > $UNMAPPED_DIR/negatives_unmapped_sort.bam
$SINGC samtools fastq $UNMAPPED_DIR/negatives_unmapped_sort.bam > $UNMAPPED_DIR/negatives_unmapped.fastq
singularity run -B $MOUNT_DIR $IMG/parasail_2.4.3.sif -t $THREADS -a sw_trace -f $UNMAPPED_DIR/negatives_unmapped.fastq -q $OLIGO_ONT_ADAPT -g $UNMAPPED_DIR/negatives_sw.sam -O {SAM}

# analyse failed reads
cat $NEGATIVES/fastq_fail/*.fastq > $WD/negatives/fastq/fail.fastq
singularity run -B $MOUNT_DIR $IMG/parasail_2.4.3.sif -t $THREADS -a sw_trace -f $WD/negatives/fastq/fail.fastq -q $OLIGO_ONT_ADAPT -g $WD/negatives/alignments_to_genome/negatives_fail.sam -O {SAM}





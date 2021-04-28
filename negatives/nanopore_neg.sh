#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

BASEDIR="/hpcnfs/scratch/FN/TL/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
singpore="/hpcnfs/scratch/TSSM/cugolini/tools/porechop/porechop.simg"
singfastp="/hpcnfs/scratch/TSSM/cugolini/tools/fastp/fastp.simg"
singtot="/hpcnfs/scratch/TSSM/cugolini/cov/img/recappable.simg"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NEGATIVES="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210423_1316_X3_FAL83184_52a41d7d/S33624_Calu3-CoV2-infected"
OLIGO_ONT_ADAPT="/hpcnfs/scratch/FN/TL/cugolini/cov/data/ont_adapter/ont_adapter.fa"

#### 	SCRIPT FOR THE ANALYSIS OF NEGATIVES SAMPLES

# pycoQC

#source activate /hpcnfs/home/ieo5215/miniconda/envs/pycoQC
#mkdir -p $WD/NEGATIVES/PYCOQC
#pycoQC -f $NEGATIVES/sequencing_summary_FAL83184_c1102ce7.txt -o $WD/NEGATIVES/PYCOQC/sequencing_summary_NEGATIVES.html
#conda deactivate

# process fastq

#mkdir -p $DATA/fastq_negatives
#cat $NEGATIVES/fastq_pass/*.fastq > $DATA/fastq_negatives/negatives.fastq
mkdir -p $WD/NEGATIVES/map_to_genome/
#singularity exec -B /hpcnfs/scratch/ $singpore minimap2 -k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7 -t 8 $GENOME_FA $DATA/fastq_negatives/negatives.fastq > $WD/NEGATIVES/map_to_genome/negatives.sam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -Sb $WD/NEGATIVES/map_to_genome/negatives.sam > $WD/NEGATIVES/map_to_genome/negatives.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -F 4 $WD/NEGATIVES/map_to_genome/negatives.bam > $WD/NEGATIVES/map_to_genome/negatives_filtered.bam 
#singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/NEGATIVES/map_to_genome/negatives_filtered.bam > $WD/NEGATIVES/map_to_genome/negatives_sorted.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/NEGATIVES/map_to_genome/negatives_sorted.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -f 16 $WD/NEGATIVES/map_to_genome/negatives.bam > $WD/NEGATIVES/map_to_genome/negatives_minus_filtered.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools sort $WD/NEGATIVES/map_to_genome/negatives_minus_filtered.bam > $WD/NEGATIVES/map_to_genome/negatives_minus_sorted.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools index $WD/NEGATIVES/map_to_genome/negatives_minus_sorted.bam

#singularity exec -B /hpcnfs/scratch/ $singpore samtools view -h -f 4 $WD/NEGATIVES/map_to_genome/negatives.bam > $WD/NEGATIVES/map_to_genome/negatives_unmapped.bam
#singularity exec -B /hpcnfs/scratch/ $singpore samtools fastq $WD/NEGATIVES/map_to_genome/negatives_unmapped.bam > $WD/NEGATIVES/map_to_genome/negatives_unmapped.fastq
singularity run -B /hpcnfs/scratch/ docker://tleonardi/parasail:2.4.3 -t 10 -a sw_trace -f $WD/NEGATIVES/map_to_genome/negatives_unmapped.fastq -q $OLIGO_ONT_ADAPT -g $WD/NEGATIVES/map_to_genome/fail/negatives_fail.sam -O {SAM}




mkdir -p $WD/NEGATIVES/map_to_genome/fail
#cat $NEGATIVES/fastq_fail/*.fastq > $DATA/fastq_negatives/negatives_fail.fastq
#singularity run -B /hpcnfs/scratch/ docker://tleonardi/parasail:2.4.3 -t 10 -a sw_trace -f $DATA/fastq_negatives/negatives_fail.fastq -q $OLIGO_ONT_ADAPT -g $WD/NEGATIVES/map_to_genome/fail/negatives_fail.sam -O {SAM}





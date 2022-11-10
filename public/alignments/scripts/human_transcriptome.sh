#!/bin/bash

set -eo pipefail
export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
REF_DATA="/hpcnfs/scratch/TSSM/cugolini/Datasets/HUMAN_REFERENCE"
GENOME_FA="$REF_DATA/Homo_sapiens.GRCh38.dna_sm.toplevel.fa"
GTF_ANNOT="$REF_DATA/Homo_sapiens.GRCh38.108.gtf"
MOUNT_DIR="/hpcnfs/scratch"
IMG="$BASEDIR/img"


# pull image
if [ ! -f "$IMG/nrceq_pipeline_latest.sif" ]; then
        cd $IMG
        singularity pull docker://cugolini/nrceq_pipeline:latest
fi


# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/nrceq_pipeline_latest.sif"

# Obtain transcriptome fasta

$SINGC bedparse gtf2bed $GTF_ANNOT | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > $REF_DATA/Homo_sapiens.GRCh38.98.bed   #from gtf(ensemble format to BED file)
$SINGC fastaFromBed -name -split -s -fi $GENOME_FA -bed $REF_DATA/Homo_sapiens.GRCh38.98.bed -fo - | perl -pe 's/>(.+)\(.\)$/>$1/' > $REF_DATA/transcriptome_fasta.fa      #get transcriptome fasta


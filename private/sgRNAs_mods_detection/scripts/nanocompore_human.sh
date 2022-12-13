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
MOUNT_DIR="/hpcnfs"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/Datasets/HUMAN_REFERENCE/transcriptome_fasta.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/Datasets/HUMAN_REFERENCE/Homo_sapiens.GRCh38.98.bed"
THREADS=3



# pull image for nanocompore
if [ ! -f "$IMG/nanocompore_v1.0.4.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/nanocompore:v1.0.4
fi
# singularity command
NANOCOMPORE="singularity exec -B $MOUNT_DIR $IMG/nanocompore_v1.0.4.sif"

# indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
BASECALLING="guppy_v601"
WD="$WD/$BASECALLING"
condition_1="WT"
condition_2="PUS7KD"
cell_line="CaCo2"

$NANOCOMPORE nanocompore sampcomp --file_list1 $WD/$cell_line/$condition_1/alignments_to_human_transcriptome/eventalign/collapse/out_eventalign_collapse.tsv --file_list2 $WD/$cell_line/$condition_2/alignments_to_human_transcriptome/eventalign/collapse/out_eventalign_collapse.tsv  --label1 $condition_1 --label2 $condition_2 --fasta $NANOCOMP_FA --outpath $WD/$cell_line/nanocompore/alignments_to_human_transcriptome/"$condition_1"_vs_"$condition_2"/ --overwrite --downsample_high_coverage 5000 --allow_warnings --min_coverage 30 --logit --nthreads $THREADS --bed $NANOCOMP_BED



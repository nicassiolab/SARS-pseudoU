#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
SCRIPTDIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/public/gRNA_mods_detection/scripts"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"
ALIGNMENTS="$BASEDIR/analysis/alignments"
WD="$BASEDIR/analysis/gRNA_mods"
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
THREADS=5
TEMP_DIR="/hpcnfs/scratch/temporary/camilla_FN/gRNA_mods"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"

BASECALLING="guppy_initial"
WD="$WD/$BASECALLING/pysamstats"

# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/nrceq_pipeline_latest.sif"

# environment
if [ ! -d "$ENVS/pysamstats/" ]; then
        conda create -n pysamstats
        source activate $ENVS/pysamstats/
	conda install -c bioconda pysamstats
        conda deactivate
fi

mkdir -p $WD
source activate $ENVS/pysamstats
pysamstats --fasta $GENOME_FA -t variation $BASEDIR/analysis/gRNA_mods/$BASECALLING/nanocompore/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam > $WD/grnas.tsv
awk 'BEGIN{OFS="\t"}{if($3=="T" && ($16/$4)>=0.05) print $1,$2,$3,$4,$16,($16/$4)}' $WD/grnas.tsv > $WD/grnas_T_C_mismatches.txt

pysamstats --fasta $GENOME_FA -t variation $BASEDIR/analysis/alignments/guppy_initial/IVT/WT/alignments_to_genome/IVT_sorted.bam > $WD/IVT.tsv
awk 'BEGIN{OFS="\t"}{if($3=="T" && ($16/$4)>=0.05) print $1,$2,$3,$4,$16,($16/$4)}' $WD/grnas.tsv > $WD/IVT_T_C_mismatches.txt



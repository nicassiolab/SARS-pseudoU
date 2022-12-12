#!/usr/bin/env bash

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov" 	
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"				# directory of configuration files used for the analysis
ALIGNMENTS="$BASEDIR/analysis/alignments"
WD="$BASEDIR/analysis/gRNA_mods"
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch"								# insert directory to mount on singularity images (needs to contain data and working directory)
GENOME_FA="$FILES/edited.fa"								# viral genome fasta
TEMP_DIR="/hpcnfs/scratch/temporary/camilla_FN/gRNA_mods"				# insert temporary directory where to copy fasta file of all samples


# parameters
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
THREADS=5
BASECALLING="guppy_initial"  #indicate basecalling version (the first basecalling used in the analysis has mixed version therefore we just call it guppy_initial)
condition="WT"


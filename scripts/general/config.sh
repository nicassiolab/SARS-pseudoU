#!/usr/bin/env bash

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/SARS-CoV-2_pseudoU/scripts/files" 	                # directory of configuration files used for the analysis
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"							# path for anaconda environments
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch"                                     	                        # insert directory to mount on singularity images (needs to contain data and working directory)
GENOME_FA="$FILES/edited.fa"                                             	                # viral genome fasta
GENOME_FAI="$FILES/edited.fa.fai"								# viral genome fasta index
TRANSCRIPTOME_ASSEMBLY="$FILES/consensus_extracted.fa"						# viral transcriptome NRCeq assembly from https://doi.org/10.1093/nar/gkac144
TRANSCRIPTOME_ASSEMBLY_FAI="$FILES/consensus_extracted.fa.fai"
NANOCOMP_BED="$FILES/aln_consensus_name_commas.bed"						# NRCeq assembly bedfile used for nanocompore processing


# parameters
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7" # mapping parameters to align viral reads to genome
TRANSCRIPTOME_PARAM="-ax map-ont -p 0 -N 10"							# mapping parameters to align reads to both viral and human transcriptome



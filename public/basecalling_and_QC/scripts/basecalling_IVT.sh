#!/bin/bash
#set -e -o pipefail
#export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

IVT="/hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Kim_2020/fast5/IVT/fast5_uncompressed"
MOUNT_DIR="/hpcnfs/scratch/"
IMG="/hpcnfs/scratch/TSSM/cugolini/cov/img"
SINGC="singularity run -B $MOUNT_DIR $IMG/guppy_latest.sif"
WD="/hpcnfs/scratch/temporary/camilla_FN/basecalling/cov"
THREADS=10
HAC_FILE="/hpcnfs/scratch/TSSM/cugolini/tools/ont-guppy-cpu/data/rna_r9.4.1_70bps_hac.cfg"
MIN_Q_SCORE=7
KIT="SQK-RNA002"
FLOWCELL="FLO-MIN106"

# basecalling with last guppy version of sample IVT 

mkdir -p $WD/IVT/fast5
cp $IVT/*.fast5 $WD/IVT/fast5
$SINGC guppy_basecaller -i $WD/IVT/fast5 -s $WD --num_callers $THREADS -c $HAC_FILE --recursive --min_qscore $MIN_Q_SCORE --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna



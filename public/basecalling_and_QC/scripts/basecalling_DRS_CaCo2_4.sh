#!/bin/bash
#set -e -o pipefail
#export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

PUS7_KD="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X2_FAL77483_43bd693a/S35756_CaCo2_C37_plus_doxy_Pus7_Pus7L_KD"
WT="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X1_FAL77093_8c3d8e3b/S35755_CaCo2_C37"
MOUNT_DIR="/hpcnfs/scratch/"
IMG="/hpcnfs/scratch/TSSM/cugolini/cov/img"
SINGC="singularity run -B $MOUNT_DIR $IMG/guppy_latest.sif"
WD="/hpcnfs/scratch/temporary/camilla_FN/basecalling/cov"
THREADS=10
HAC_FILE="/hpcnfs/scratch/TSSM/cugolini/tools/ont-guppy-cpu/data/rna_r9.4.1_70bps_hac.cfg"
MIN_Q_SCORE=7
KIT="SQK-RNA002"
FLOWCELL="FLO-MIN106"

# basecalling with last guppy version of sample DRS_CaCo2_4, WT and KD 

mkdir -p $WD/c37_WT/fast5
cp $WT/fast5_pass/*.fast5 $WD/c37_WT/fast5
cp $WT/fast5_fail/*.fast5 $WD/c37_WT/fast5
$SINGC guppy_basecaller -i $WD/c37_WT/fast5 -s $WD --num_callers $THREADS -c $HAC_FILE --recursive --min_qscore $MIN_Q_SCORE --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna

mkdir -p $WD/c37_PUS7KD/fast5
cp $PUS7_KD/fast5_pass/*.fast5 $WD/c37_PUS7KD/fast5
cp $PUS7_KD/fast5_fail/*.fast5 $WD/c37_PUS7KD/fast5

$SINGC guppy_basecaller -i $WD/c37_PUS7KD/fast5 -s $WD --num_callers $THREADS -c $HAC_FILE --recursive --min_qscore $MIN_Q_SCORE --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna


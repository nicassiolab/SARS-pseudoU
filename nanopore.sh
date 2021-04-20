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
TRANSCRIPTOME_ASSEMBLY="$BASEDIR/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_transcripts.fas"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
PUS7_KD="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210414_1407_X2_FAL77079_62af908f/S33250_CaCo2_C34_plus_doxy_Pus7KD"
WT="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210414_1407_X1_FAL82866_fb16777a/S33249_CaCo2_C34"

#### 	SCRIPT FOR THE ANALYSIS OF PUS7 KD SAMPLES

mkdir $DATA/fastq_PUS7
cat $PUS7_KD/fastq_pass/*.fastq > $DATA/fastq_PUS7/PUS7_KD.fastq
cat $WT/fastq_pass/*.fastq > $DATA/fastq_PUS7/WT.fastq

/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux index -d $PUS7_KD/fast5_pass $DATA/fastq_PUS7/PUS7_KD.fastq
/hpcnfs/scratch/TSSM/cugolini/tools/f5c/f5c-v0.6/f5c_x86_64_linux index -d $WT/fast5_pass $DATA/fastq_PUS7/WT.fastq

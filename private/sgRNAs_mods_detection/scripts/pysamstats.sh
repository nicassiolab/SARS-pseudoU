#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
SCRIPTDIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/public/gRNA_mods_detection/scripts"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/files"
ALIGNMENTS="$BASEDIR/analysis/alignments"
WD="$BASEDIR/analysis/sgRNAs_mods_detection"
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_PARAM="-k 8 -w 1 -ax splice -g 30000 -G 30000 -A1 -B2 -O2,24 -E1,0 -C0 -z 400,200 --no-end-flt -F 40000 -N 32 --splice-flank=no --max-chain-skip=40 -un -p 0.7"
THREADS=5
TEMP_DIR="/hpcnfs/scratch/temporary/camilla_FN/gRNA_mods"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"

BASECALLING="guppy_v601"
CELL_LINE="CaCo2"
WD="$WD/$BASECALLING/$CELL_LINE/pysamstats"

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
for condition in WT PUS7KD; do 
	pysamstats --fasta $TRANSCRIPTOME_ASSEMBLY -t variation $ALIGNMENTS/$BASECALLING/$CELL_LINE/$condition/alignments_to_assembly/DRS_CaCo2_4_nanocompore.bam > $WD/"$condition".tsv
	awk 'BEGIN{OFS="\t"}{if($3=="T" && ($16/$4)>=0.05) print $1,$2,$3,$4,$16,($16/$4)}' $WD/"$condition".tsv > $WD/"$condition"_T_C_mismatches.txt
done




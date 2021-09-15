#!/bin/bash
#set -e -o pipefail
#export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

BASEDIR="/hpcnfs/scratch/FN/TL/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
singpore="/hpcnfs/scratch/TSSM/cugolini/tools/porechop/porechop.simg"
singfastp="/hpcnfs/scratch/TSSM/cugolini/tools/fastp/fastp.simg"
singtot="/hpcnfs/scratch/TSSM/cugolini/cov/img/recappable.simg"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
PUS7_KD="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210414_1407_X2_FAL77079_62af908f/S33250_CaCo2_C34_plus_doxy_Pus7KD"
WT="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210414_1407_X1_FAL82866_fb16777a/S33249_CaCo2_C34"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/consensus_extracted.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"



source activate $ENVS/nanocompore_v100rc32/

for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/eventalign/IVT/eventalign*; do
        base=${filename##*/eventalign_}
	/hpcnfs/home/ieo5215/miniconda/envs/nanocompore_v100rc32/bin/nanocompore sampcomp --file_list1 $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/eventalign/WT/out_eventalign_collapse.tsv --file_list2 "$filename"/out_eventalign_collapse.tsv  --label1 WT --label2 IVT --fasta $NANOCOMP_FA --outpath $WD/PUS7_KD/map_to_recap_assembly/NANOCOMPORE/sampcomp/WT/"$base" --overwrite --downsample_high_coverage 5000 --allow_warnings --pvalue_thr 0.01 --min_coverage 30 --logit --nthreads 3 --bed $NANOCOMP_BED
done



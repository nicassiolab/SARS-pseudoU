#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

BASEDIR="/hpcnfs/scratch/FN/TL/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
NANORMS="/hpcnfs/scratch/TSSM/cugolini/tools/nanoRMS"
singpore="/hpcnfs/scratch/TSSM/cugolini/tools/porechop/porechop.simg"
singfastp="/hpcnfs/scratch/TSSM/cugolini/tools/fastp/fastp.simg"
singtot="/hpcnfs/scratch/TSSM/cugolini/cov/img/recappable.simg"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
PUS7_KD="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X2_FAL77483_43bd693a/S35756_CaCo2_C37_plus_doxy_Pus7_Pus7L_KD"
WT="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X1_FAL77093_8c3d8e3b/S35755_CaCo2_C37"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/consensus_extracted.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"

mkdir -p $WD/PUS7_KD_C37/nanoRMS
source activate $ENVS/nanoRMS/
#cp $TRANSCRIPTOME_ASSEMBLY $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta
#java -jar $NANORMS/epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=$WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta OUTPUT=$WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta.dict
#samtools faidx $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta
#cd $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/
#python3 $NANORMS/epinano_RMS/epinano_rms.py -R $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta -b WT_sorted.bam -s $NANORMS/epinano_RMS/sam2tsv.jar -n 4
#python3 $NANORMS/epinano_RMS/epinano_rms.py -R $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta -b PUS7_KD_sorted.bam -s $NANORMS/epinano_RMS/sam2tsv.jar -n 4
#conda deactivate
#source activate $ENVS/R_env/
#cd $NANORMS/predict_rna_mod
#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_pairedcondition_transcript.R -f $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/WT_sorted.per.site.baseFreq.csv -s $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/PUS7_KD_sorted.per.site.baseFreq.csv 

#### calculate stoichiometry
cd $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/
## Convert BAM to TSV
#samtools view -h WT_sorted.bam | java -jar $NANORMS/epinano_RMS/sam2tsv.jar -r $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta  > WT.bam.tsv 
#samtools view -h PUS7_KD_sorted.bam | java -jar $NANORMS/epinano_RMS/sam2tsv.jar -r $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta  > PUS7_KD.bam.tsv
## Run EpiNano
python3.6 TSV_to_Variants_Freq.py3 -f WT.bam.tsv -t 10 


#!/bin/bash
set -e -o pipefail

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
ANALYSIS="$BASEDIR/analysis"
SCRIPTDIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/public/gRNA_mods_detection/scripts"
FILES="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/public/gRNA_mods_detection/files"
WD="$BASEDIR/analysis/gRNA_mods"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
NANORMS="/hpcnfs/scratch/TSSM/cugolini/tools/nanoRMS"

THREADS=5

# environment
if [ ! -d "$ENVS/nanoRMS/" ]; then
	conda create -n nanoRMS --file $FILES/nanoRMS_requirements.txt
	source activate $ENVS/nanoRMS/
	$ENVS/nanoRMS/bin/pip install sklearn=0.23.1
	conda deactivate
fi

BASECALLING="guppy_initial"
WD="$WD/$BASECALLING"


# create dictionary for viral genome
#source activate $ENVS/nanoRMS/
#mkdir -p $WD/nanoRMS
#cp $GENOME_FA $WD/nanoRMS/
#java -jar $NANORMS/epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=$WD/nanoRMS/edited.fa OUTPUT=$WD/nanoRMS/edited.fa.dict
#samtools faidx $WD/nanoRMS/edited.fa
#conda deactivate


# predict RNA mods 
source activate $ENVS/nanoRMS/
cd $WD/nanoRMS

#cp $WD/nanocompore/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam* $NANORMS/epinano_RMS/
#python3 $NANORMS/epinano_RMS/epinano_rms.py -R $WD/guppy_initial/nanoRMS/edited.fa -b gRNAs_sorted.bam -s $NANORMS/epinano_RMS/sam2tsv.jar -n $THREADS
#mv $NANORMS/epinano_RMS/gRNAs_sorted* $WD/nanoRMS/

cp $ANALYSIS/alignments/guppy_initial/IVT/WT/alignments_to_genome/IVT_sorted.bam* $WD/nanoRMS
#python3 $NANORMS/epinano_RMS/epinano_rms.py -R $WD/guppy_initial/nanoRMS/edited.fa -b IVT_sorted.bam -s $NANORMS/epinano_RMS/sam2tsv.jar -n $THREADS
#mv $NANORMS/epinano_RMS/IVT_sorted* $WD/guppy_initial/nanoRMS/

conda deactivate

#source activate $ENVS/R_env/
#mkdir -p $WD/nanoRMS/IVT $WD/nanoRMS/gRNAs $WD/nanoRMS/IVT_gRNAs
#cd $NANORMS/predict_rna_mod
#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_singlecondition.R -f $WD/nanoRMS/IVT_sorted.per.site.baseFreq.csv
#mv $NANORMS/predict_rna_mod/Sample1_predicted_y_sites.tsv $WD/nanoRMS/IVT
#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_singlecondition.R -f $WD/nanoRMS/gRNAs_sorted.per.site.baseFreq.csv
#mv $NANORMS/predict_rna_mod/Sample1_predicted_y_sites.tsv $WD/nanoRMS/gRNAs
#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_pairedcondition_transcript.R -f $WD/nanoRMS/IVT_sorted.per.site.baseFreq.csv -s $WD/nanoRMS/gRNAs_sorted.per.site.baseFreq.csv
#mkdir -p $WD/nanoRMS/IVT_gRNAs
#mv $NANORMS/predict_rna_mod/Paired_comparison_altering_sites_pU_predictions.* $WD/nanoRMS/IVT_gRNAs





#### calculate stoichiometry
#cd $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/
## Convert BAM to TSV
#samtools view -h WT_sorted.bam | java -jar $NANORMS/epinano_RMS/sam2tsv.jar -r $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta  > WT.bam.tsv 
#samtools view -h PUS7_KD_sorted.bam | java -jar $NANORMS/epinano_RMS/sam2tsv.jar -r $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta  > PUS7_KD.bam.tsv
## Run EpiNano
#python3.6 TSV_to_Variants_Freq.py3 -f WT.bam.tsv -t 10 








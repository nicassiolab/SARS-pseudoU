#!/bin/bash
set -e -o pipefail

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
ANALYSIS="$BASEDIR/analysis"
SCRIPTDIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/public/sgRNAs_mods_detection/scripts"
WD="$BASEDIR/analysis/sgRNAs_mods_detection"
GENOME_VIRAL="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
TRANSCRIPTOME_ASSEMBLY="$ANALYSIS/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
GENOME_HG="/hpcnfs/scratch/TSSM/cugolini/Datasets/HUMAN_REFERENCE/Homo_sapiens.GRCh38.dna_sm.toplevel.fa"
TRANSCRIPTOME_HG="/hpcnfs/scratch/TSSM/cugolini/Datasets/HUMAN_REFERENCE/transcriptome_fasta.fa"
IMG="$BASEDIR/img"
MOUNT_DIR="/hpcnfs/scratch"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
NANORMS="/hpcnfs/scratch/TSSM/cugolini/tools/nanoRMS"

THREADS=6
BASECALLING="guppy_v601"
CELL_LINE="CaCo2"
SAMPLE="DRS_CaCo2_4"
WD=$WD/$BASECALLING/$CELL_LINE

# environment
if [ ! -d "$ENVS/nanoRMS/" ]; then
	conda create -n nanoRMS --file $FILES/nanoRMS_requirements.txt
	source activate $ENVS/nanoRMS/
	$ENVS/nanoRMS/bin/pip install sklearn=0.23.1
	conda deactivate
fi

# create dictionary for viral transcriptome

#source activate $ENVS/nanoRMS/
#mkdir -p $WD/nanoRMS/viral/assembly
#cp $TRANSCRIPTOME_ASSEMBLY $WD/nanoRMS/viral/assembly
#java -jar $NANORMS/epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=$WD/nanoRMS/viral/assembly/consensus_extracted.fa OUTPUT=$WD/nanoRMS/viral/assembly/consensus_extracted.fa.dict
#samtools faidx $WD/nanoRMS/viral/assembly/consensus_extracted.fa
#conda deactivate

# create dictionary for viral genome

#source activate $ENVS/nanoRMS/
#mkdir -p $WD/nanoRMS/viral/genome
#cp $GENOME_VIRAL $WD/nanoRMS/viral/genome
#java -jar $NANORMS/epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=$WD/nanoRMS/viral/genome/edited.fa OUTPUT=$WD/nanoRMS/viral/genome/edited.fa.dict
#samtools faidx $WD/nanoRMS/viral/genome/edited.fa
#conda deactivate

# create dictionary for human transcriptome

#source activate $ENVS/nanoRMS/
#mkdir -p $WD/nanoRMS/human/transcriptome
#cp $TRANSCRIPTOME_HG $WD/nanoRMS/human/transcriptome
#java -jar $NANORMS/epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=$WD/nanoRMS/human/transcriptome/transcriptome_fasta.fa OUTPUT=$WD/nanoRMS/human/transcriptome/transcriptome_fasta.fa.dict
#samtools faidx $WD/nanoRMS/human/transcriptome/transcriptome_fasta.fa
#conda deactivate

#create dictionary for human genome

#source activate $ENVS/nanoRMS/
#mkdir -p $WD/nanoRMS/human/genome
#cp $GENOME_HG $WD/nanoRMS/human/genome
#java -jar $NANORMS/epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=$WD/nanoRMS/human/genome/Homo_sapiens.GRCh38.dna_sm.toplevel.fa OUTPUT=$WD/nanoRMS/human/genome/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.dict
#samtools faidx $WD/nanoRMS/human/genome/Homo_sapiens.GRCh38.dna_sm.toplevel.fa
#conda deactivate 

#predict RNA mods 


for organism in human ; do
	if [ $organism = "human" ]; then
		tr_ref="transcriptome"
        elif [ $organism = "viral" ]; then
		tr_ref="assembly"
	fi

	for reference in $tr_ref ; do
		if [ $reference != "genome" ]; then
			suffix="nanocompore"
		elif [ $reference = "genome" ]; then
			suffix="sorted"
		fi

		for condition in WT  ; do

#			cp $ANALYSIS/alignments/$BASECALLING/$CELL_LINE/$condition/alignments_to_"$organism"_"$reference"/"$SAMPLE"_"$suffix".bam* $NANORMS/epinano_RMS/
			cd $NANORMS/epinano_RMS/
			source activate $ENVS/nanoRMS/
			python3 $NANORMS/epinano_RMS/epinano_rms.py -R $WD/nanoRMS/$organism/$reference/transcriptome_fasta.fa -b "$SAMPLE"_"$suffix".bam -s $NANORMS/epinano_RMS/sam2tsv.jar -n $THREADS
#			conda deactivate
#			mkdir -p $WD/nanoRMS/$condition
#			mv $NANORMS/epinano_RMS/"$SAMPLE"_sorted* $WD/nanoRMS/$condition
#			source activate $ENVS/R_env/
#			cd $NANORMS/predict_rna_mod
#			Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_singlecondition.R -f $WD/nanoRMS/$condition/"$SAMPLE"_nanocompore.per.site.baseFreq.csv
#			mv $NANORMS/predict_rna_mod/Sample1_predicted_y_sites.tsv $WD/nanoRMS/$condition/
#			conda deactivate
		done
	done
done


# perform paired conditions analysis

#source activate $ENVS/R_env/
#cd $NANORMS/predict_rna_mod/
#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_pairedcondition_transcript.R -f $WD/nanoRMS/WT/"$SAMPLE"_nanocompore.per.site.baseFreq.csv -s $WD/nanoRMS/PUS7KD/"$SAMPLE"_nanocompore.per.site.baseFreq.csv
#mv Paired* $WD/nanoRMS/
#conda deactivate




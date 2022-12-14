#!/bin/bash
set -e -o pipefail
export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
NANORMS="/hpcnfs/scratch/TSSM/cugolini/tools/nanoRMS"

singpore="/hpcnfs/scratch/TSSM/cugolini/tools/porechop/porechop.simg"
singfastp="/hpcnfs/scratch/TSSM/cugolini/tools/fastp/fastp.simg"
singtot="/hpcnfs/scratch/TSSM/cugolini/cov/img/recappable.simg"
PUS7_KD="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X2_FAL77483_43bd693a/S35756_CaCo2_C37_plus_doxy_Pus7_Pus7L_KD"
WT="/hpcnfs/techunits/genomics/PublicData/TSSM/tleonardi/FAST5/20210825_1506_X1_FAL77093_8c3d8e3b/S35755_CaCo2_C37"
TRANSCRIPTOME_ASSEMBLY="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/recappable_assembly/two_datasets/assemblies/pinfish/consensus_extraction/consensus_extracted.fa"
NANOCOMP_FA="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/consensus_extracted.fa"
NANOCOMP_BED="/hpcnfs/scratch/TSSM/cugolini/cov/scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus_name_commas.bed"
IVT="/hpcnfs/scratch/TSSM/tleonardi/SARS-CoV-2-datasets/Kim_2020/fastq/IVT/ont_research.fa"

#	script to run nanoRMS on SARS-CoV-2 data from sample C37

mkdir -p $WD/PUS7_KD_C37/nanoRMS
source activate $ENVS/nanoRMS/
#cp $TRANSCRIPTOME_ASSEMBLY $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta
#java -jar $NANORMS/epinano_RMS/picard.jar CreateSequenceDictionary REFERENCE=$WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta OUTPUT=$WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta.dict
#samtools faidx $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta
#conda deactivate

######### PREDICT RNA MODS (SINGLE FILE) ON WT and PUS7KD (SINGLE FILE) MAPPED TO TRANSCRIPTOME


###	IVT
mkdir -p /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/IVT
mkdir -p /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT_vs_IVT

for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/map_to_transcriptome/IVT/*.bam; do
	name=${filename##*/}
	base=${name%_sorted.bam}
	if [ ! -e /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/IVT/"$base"_sorted.per.site.baseFreq.csv ]; then
		cp $filename $NANORMS/epinano_RMS/
		source activate $ENVS/nanoRMS/
		cd $NANORMS/epinano_RMS/
		python3 $NANORMS/epinano_RMS/epinano_rms.py -R $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta -b "$base"_sorted.bam -s $NANORMS/epinano_RMS/sam2tsv.jar -n 12
		cp -r $NANORMS/epinano_RMS/"$base"_sorted.tmp_splitted_base_freq/ /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/IVT/"$base"_sorted.tmp_splitted_base_freq/
		mv $NANORMS/epinano_RMS/"$base"_sorted.per.site.baseFreq.csv /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/IVT/"$base"_sorted.per.site.baseFreq.csv
		rm $NANORMS/epinano_RMS/"$base"_sorted.bam*
		rm -r $NANORMS/epinano_RMS/"$base"_sorted.tmp_splitted_base_freq/
		conda deactivate
	fi
	#source activate $ENVS/R_env/
	#cd $NANORMS/predict_rna_mod
	#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_singlecondition.R -f /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/IVT/"$base"_sorted.per.site.baseFreq.csv
	#mv $NANORMS/predict_rna_mod/Sample1_predicted_y_sites.tsv /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/IVT
	#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_pairedcondition_transcript.R -f /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT/WT_sorted.per.site.baseFreq.csv -s /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/IVT/"$base"_sorted.per.site.baseFreq.csv
	#mkdir -p /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT_vs_IVT/"$base"
	#mv $NANORMS/predict_rna_mod/Paired_comparison_altering_sites_pU_predictions.* /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT_vs_IVT/"$base"
	#conda deactivate	
done


#source activate $ENVS/nanoRMS/
#cd $NANORMS/epinano_RMS/

###	WT
#mv /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT_sorted.bam* $NANORMS/epinano_RMS/
#python3 $NANORMS/epinano_RMS/epinano_rms.py -R $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta -b WT_sorted.bam -s $NANORMS/epinano_RMS/sam2tsv.jar -n 12
#mv $NANORMS/epinano_RMS/WT_sorted* /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/


###	PUS7_KD
#mv /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/PUS7_KD/PUS7_KD_sorted.bam* $NANORMS/epinano_RMS/
#python3 $NANORMS/epinano_RMS/epinano_rms.py -R $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta -b PUS7_KD_sorted.bam -s $NANORMS/epinano_RMS/sam2tsv.jar -n 12
#mv $NANORMS/epinano_RMS/PUS7_KD_sorted.bam* /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/PUS7_KD
#cp -r $NANORMS/epinano_RMS/PUS7_KD_sorted.tmp_splitted_base_freq/ /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/PUS7_KD
#rm -r $NANORMS/epinano_RMS/PUS7_KD_sorted.tmp_splitted_base_freq/
#mv $NANORMS/epinano_RMS/PUS7_KD_sorted.per.site.baseFreq.csv /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/PUS7_KD


#conda deactivate


#source activate $ENVS/R_env/
#cd $NANORMS/predict_rna_mod
#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_singlecondition.R -f /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT_sorted.per.site.baseFreq.csv
#mv $NANORMS/predict_rna_mod/Sample1_predicted_y_sites.tsv /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT
#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_singlecondition.R -f /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/PUS7_KD/PUS7_KD_sorted.per.site.baseFreq.csv
#mv $NANORMS/predict_rna_mod/Sample1_predicted_y_sites.tsv /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/PUS7_KD
#Rscript --vanilla $NANORMS/predict_rna_mod/Pseudou_prediction_pairedcondition_transcript.R -f /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT/WT_sorted.per.site.baseFreq.csv -s /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/PUS7_KD/PUS7_KD_sorted.per.site.baseFreq.csv 
#mkdir -p /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT_vs_PUS7_KD
#mv $NANORMS/predict_rna_mod/Paired_comparison_altering_sites_pU_predictions.* /hpcnfs/scratch/temporary/camilla_TL/temp_dir_nanop/WT_vs_PUS7_KD



#### calculate stoichiometry
#cd $WD/PUS7_KD_C37/map_to_recap_assembly/NANOCOMPORE/map_to_transcriptome/
## Convert BAM to TSV
#samtools view -h WT_sorted.bam | java -jar $NANORMS/epinano_RMS/sam2tsv.jar -r $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta  > WT.bam.tsv 
#samtools view -h PUS7_KD_sorted.bam | java -jar $NANORMS/epinano_RMS/sam2tsv.jar -r $WD/PUS7_KD_C37/nanoRMS/consensus_extracted.fasta  > PUS7_KD.bam.tsv
## Run EpiNano
#python3.6 TSV_to_Variants_Freq.py3 -f WT.bam.tsv -t 10 








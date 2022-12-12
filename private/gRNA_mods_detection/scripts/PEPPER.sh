#!/bin/bash

BASE="/hpcnfs/scratch/TSSM/cugolini/cov/analysis"
GENOME_FA="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa"
GENOME_FAI="/hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa.fai"
BASECALLING="guppy_v601"
CELL_LINE="CaCo2"
WD="$BASE/sgRNAs_mods_detection/BASECALLING/CELL_LINE/PEPPER"


# Set up input data
INPUT_DIR="$WD"
SAMPLE="DRS_CaCo2_4"
REF="edited.fa"


# Set the number of CPUs to use
THREADS="8"

for condition in WT PUS7KD; do

	# Set up output directory
	OUTPUT_DIR="$WD/output/$condition"
	OUTPUT_PREFIX="SARS_CoV_2_sgRNAs_GRCh38_PEPPER_Margin_DeepVariant"

	## Create local directory structure
	mkdir -p "${OUTPUT_DIR}"
	mkdir -p "${INPUT_DIR}"

	# Download the data to input directory
	cp ${BASE}/alignments/$BASECALLING/$CELL_LINE/$condition/alignments_to_genome/"$SAMPLE"_sorted.bam ${INPUT_DIR}
	cp ${BASE}/alignments/$BASECALLING/$CELL_LINE/$condition/alignments_to_genome/"$SAMPLE"_sorted.bam.bai ${INPUT_DIR}
	cp $GENOME_FA ${INPUT_DIR}
	cp $GENOME_FAI ${INPUT_DIR}

done

# Pull the docker images
singularity pull docker://jmcdani20/hap.py:v0.3.12
singularity pull docker://kishwars/pepper_deepvariant:r0.8

# Run PEPPER-Margin-DeepVariant
#singularity exec --bind /usr/lib/locale/ \
#pepper_deepvariant_r0.8.sif \
#run_pepper_margin_deepvariant call_variant \
#-b "${INPUT_DIR}/${BAM}" \
#-f "${INPUT_DIR}/${REF}" \
#-o "${OUTPUT_DIR}" \
#-p "${OUTPUT_PREFIX}" \
#-t ${THREADS} \
#-r NC_045512v2:1-29003 


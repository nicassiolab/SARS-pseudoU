BASE="/hpcnfs/scratch/TSSM/cugolini/cov/analysis/gRNA_mods"

# Set up input data
INPUT_DIR="${BASE}/all_cell_lines/PEPPER/data"
REF="edited.fa"
BAM="gRNAs_sorted.bam"

# Set the number of CPUs to use
THREADS="8"

# Set up output directory
OUTPUT_DIR="${BASE}/all_cell_lines/PEPPER/output"
OUTPUT_PREFIX="SARS_CoV_2_gRNAs_GRCh38_PEPPER_Margin_DeepVariant"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"

# Download the data to input directory
cp ${BASE}/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam ${INPUT_DIR}
cp ${BASE}/all_cell_lines/alignments_to_genome/gRNAs_sorted.bam.bai ${INPUT_DIR}
cp /hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa ${INPUT_DIR}
cp /hpcnfs/scratch/TSSM/cugolini/CoV-2_analysis/reference_genome/results/edited.fa.fai ${INPUT_DIR}

# Pull the docker images
singularity pull docker://jmcdani20/hap.py:v0.3.12
singularity pull docker://kishwars/pepper_deepvariant:r0.8

# Run PEPPER-Margin-DeepVariant
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.8.sif \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t ${THREADS} \
-r NC_045512v2:1-29003 


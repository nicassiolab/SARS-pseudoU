library(tidyverse)
library(dplyr)
library(ggpubr)
library(tidyr)
library(R.utils)
library(reshape2)
library(rlang)
library(xlsx)
library(tidyr)
library(stringr)
library(xlsx)
library(seqinr)
library(DescTools)
library(Biostrings)



# Directories to set
ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"                                    # insert directory where the analysis will be stored
CURRDIR="/Volumes/scratch/FN/TL/cugolini/cov/SARS-CoV-2_pseudoU/scripts"        # insert script directory
which_nucl <- "U"                                                               # nucleotide on which to look for the modification 
BASECALLING="guppy_initial"                                                     # guppy version used for basecalling

#Source functions and set directories 
source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))     # source script with functions
FILES=cdp("files")                                                              # directories where useful files for the analysis are found
RESULTSDIR=bdp("results/sgRNAs_mods_detection_new_coords")                      # directory where results will be stores
IVT_bedfile=fdp("IVT_junctions.bed")                                            # IVT bedfile from Kim et al.
assembly_bedfile=fdp("aln_consensus_name_commas.bed")                           # NRCeq assembly bedfile from Ugolini et al.
ORF_annotation_bedfile=fdp("orf_annotate.bed")                                  # ORF annotation from Ugolini et al.
fragments_bedfile=fdp("fragments_1_based.txt")                                  # genomic fragments developed for the analysis
viral_genome_fasta=fdp("edited.fa")                                             # reference viral genome fasta edited as in Kim et al.
pseudoU_motifs=fdp("sites_pseudoU_motifs_blacklist.txt")                        # list of putative motifs associate to pseudouridylation

#Datasets
cell_line <- as.list(c("CaLu3","CaCo2","VeroE6"))                               # list of cell lines included in the analysis

calu3_results <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/CaLu3/nanocompore/sampcomp/"))
caco2_results <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/CaCo2/nanocompore/sampcomp/"))
vero_results <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/VeroE6/nanocompore/sampcomp/"))

all_cell_lines_results <- "/Volumes/scratch/temporary/cugolini/nanocompore/sampcomp/" # directory with nanocompore results from all cell lines pulled together

########################### PARAMETERS #########################################

LOR_thresh <- 0.5                                                               # Log Odds Ratio threshold                                                          
pval_thresh <- 0.01                                                             # GMM Logit p-value threshold
IVT_junc_interval_left <- 25                                                    # numbers of nucleotides upstream of an IVT junction point to include as part of the junction 
IVT_junc_interval_right <- 25                                                   # numbers of nucleotides downstream of an IVT junction point to include as part of the junction 
ORF_junc_interval_left <- 15                                                    # numbers of nucleotides upstream of an ORF junction point to include as part of the junction 
ORF_junc_interval_right <- 15                                                   # numbers of nucleotides downstream of an ORF junction point to include as part of the junction 
fleming_paper_notation <- c(22322,23317,27164,28417,28759,28927,29418)          # sites called as pseudouridines in Fleming et al. paper
fleming_sites_bed <- fleming_paper_notation-3                                   # move Fleming sites to Nanocompore notation

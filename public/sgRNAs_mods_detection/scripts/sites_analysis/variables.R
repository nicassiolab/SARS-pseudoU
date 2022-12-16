library(tidyverse)
library(tidyr)
library(R.utils)
library(reshape2)
library(rlang)
library(stringr)
library(dplyr)
library(xlsx)
library(seqinr)
library(ggpubr)
library(DescTools)


# Directories to set
ROOTDIR=""                                                                      # insert directory where the analysis will be stored
CURRDIR=""                                                                      # insert script directory
which_nucl <- "U"                                                               # nucleotide on which to look for the modification 
BASECALLING="guppy_initial"                                                     # guppy version used for basecalling

#Source functions and set directories 
source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))     # source script with functions
FILES=cdp("files")                                                              # directories where useful files for the analysis are found
RESULTSDIR=bdp("results/sgRNAs_mods_detection")                                 # directory where results will be stores
IVT_bedfile=fdp("IVT_junctions.bed")                                            # IVT bedfile from Kim et al.
assembly_bedfile=fdp("aln_consensus_name_commas.bed")                           # NRCeq assembly bedfile from Ugolini et al.
ORF_annotation_bedfile=fdp("orf_annotate.bed")                                  # ORF annotation from Ugolini et al.
fragments_bedfile=fdp("fragments_genomic_coord_UCSC.txt")                       # genomic fragments developed for the analysis
viral_genome_fasta=fdp("edited.fa")                                             # reference viral genome fasta edited as in Kim et al.
pseudoU_motifs=fdp("sites_pseudoU_motifs_blacklist.txt")                        # list of putative motifs associate to pseudouridylation

#Datasets
cell_line <- as.list(c("calu3","caco2","vero")) # list of cell lines included in the analysis
calu3_results <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/CaLu3/nanocompore/sampcomp/"))
caco2_results <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/CaCo2/nanocompore/sampcomp/"))
vero_results <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/VeroE6/nanocompore/sampcomp/"))

########################### PARAMETERS #########################################

LOR_thresh <- 0.5                                                               # Log Odds Ratio threshold                                                          
pval_thresh <- 0.01                                                             # GMM Logit p-value threshold
IVT_junc_interval_left <- 25                                                    # numbers of nucleotides upstream of an IVT junction point to include as part of the junction 
IVT_junc_interval_right <- 25                                                   # numbers of nucleotides downstream of an IVT junction point to include as part of the junction 
ORF_junc_interval_left <- 15                                                    # numbers of nucleotides upstream of an ORF junction point to include as part of the junction 
ORF_junc_interval_right <- 15                                                   # numbers of nucleotides downstream of an ORF junction point to include as part of the junction 
fleming_paper_notation <- c(22322,23317,27164,28417,28759,28927,29418)          # sites called as pseudouridines in Fleming et al. paper
fleming_sites_bed <- fleming_paper_notation-3                                   # move Fleming sites to Nanocompore notation

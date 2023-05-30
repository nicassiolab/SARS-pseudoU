
ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
CURRDIR="/Volumes/scratch/FN/TL/cugolini/cov/SARS-CoV-2_pseudoU/scripts"                                                                      # insert script directory
which_nucl <- "U"                                                               # nucleotide on which to look for the modification 
BASECALLING="guppy_initial"                                                     # guppy version used for basecalling

# source functions 
source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))

RESULTSDIR=bdp("results/gRNAs_mods_detection_2100")
nanocompore_gRNAs_file <- bdp(paste0("analysis/gRNA_mods/",
                                     BASECALLING,
                                     "/nanocompore/all_cell_lines/nanocompore_downsample_2100/outnanocompore_results.tsv"))
SGRNASDIR <- bdp("results/sgRNAs_mods_detection_new_coords/shared/bedtracks")
SGRNA_FINAL <- bdp("results/sgRNAs_mods_detection_new_coords/shared_U.xls")
SGRNAS_all_cell_lines <- bdp(paste0("results/sgRNAs_mods_detection_new_coords/","all_cell_lines_modified_sites",which_nucl,".xls"))
FILES=cdp("files")                                                              # directories where useful files for the analysis are found

IVT_bedfile=fdp("IVT_junctions.bed")                                            # IVT bedfile from Kim et al.
fragments_bedfile=fdp("fragments_1_based.txt")                                  # genomic fragments developed for the analysis
SNPs_file <- fdp("SNPs.txt")
pseudoU_motifs=fdp("sites_pseudoU_motifs_blacklist.txt")                        # list of putative motifs associate to pseudouridylation

#Pysamstats results
calu3_pysam <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/CaLu3/pysamstats/alignments_to_genome/"))
caco2_pysam <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/CaCo2/pysamstats/alignments_to_genome/"))
vero_pysam <- bdp(paste0("analysis/sgRNAs_mods_detection/",BASECALLING,"/VeroE6/pysamstats/alignments_to_genome/"))


########################### PARAMETERS #########################################

LOR_thresh <- 0.5                                                               # insert LOR threshold
pval_thresh <- 0.01                                                             # insert GMM logit pvalue threshold
pval_thresh_gRNAs <- 0.01                                                       # insert GMM logit pvalue threshold for the gRNAs analysis
IVT_junc_interval_left <- 25                                                    # insert number of nucleotide from the left extreme of an IVT junction to consider part of the junction 
IVT_junc_interval_right <- 25                                                   # insert number of nucleotide from the right extreme of an IVT junction to consider part of the junction 
fleming_paper_notation <- c(22322,23317,27164,28417,28759,28927,29418)          # sites present in reference paper Fleming et al.
fleming_sites_bed <- fleming_paper_notation-3                                   # move Fleming et al. notation to Nanocompore notation
SPIKETRSB <- 21555
nucl_n_from_spike_TRSB_to_end <- (29902-SPIKETRSB)                              # genome length after SPIKE TRS-B used for hypergeometric test calculation
kmer_size <- 5                                                                  # size of k-mers in nanocompore

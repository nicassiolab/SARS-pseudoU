source(paste(dirname(getSourceEditorContext()$path),"functions.R",sep="/"))

ROOTDIR="/Volumes/scratch/TSSM/cugolini/cov"
RESULTSDIR=bdp("analysis/gRNA_mods/guppy_initial/nanocompore/downstream")
nanocompore_gRNAs_file <- bdp("analysis/gRNA_mods/guppy_initial/nanocompore/all_cell_lines/nanocompore/outnanocompore_results.tsv")
SGRNASDIR <- ("/Volumes/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/shared/bedtracks")
IVT_junctions_file <- bdp("scripts_new/backupped_data/IVT_junctions.bed")
assembly_file <- bdp("scripts_new/backupped_data/data_huxelerate_extraction/transcriptome_assembly/aln_consensus.bed")
fragments_file <- bdp("scripts_new/backupped_data/RAPID/fragments_genomic_coord_UCSC.txt")
SNPs_file <- ("/Volumes/scratch/FN/TL/cugolini/cov/scripts/files/SNPs.txt")
sitelist_pseudoU_file <- bdp("scripts_new/backupped_data/sites_pseudoU_motifs_blacklist.txt")
sitelist_m6A_file <- bdp("scripts_new/backupped_data/sites_blacklist.txt")

########################### PARAMETERS #########################################

LOR_thresh <- 0.5                                                               # insert LOR threshold
pval_thresh <- 0.01                                                             # insert GMM logit pvalue threshold
IVT_junc_interval_left <- 25                                                    # insert number of nucleotide from the left extreme of an IVT junction to consider part of the junction 
IVT_junc_interval_right <- 25                                                   # insert number of nucleotide from the right extreme of an IVT junction to consider part of the junction 
fleming_paper_notation <- c(22322,23317,27164,28417,28759,28927,29418)          # sites present in reference paper Fleming et al.
fleming_sites_bed <- fleming_paper_notation-3                                   # move Fleming et al. notation to Nanocompore notation
which_nucl <- "U"                                                               # nucleotide for which detecting mods: possibilities are U or A


Description of the directory

This directory contains analysis for modifications in three different cell lines infected with SARS-CoV-2 : 
CaCo2 - 4 datasets : WTC37,WTC34,SRAFFAELE(paper assembly),MATTHEWS
CaLu3 - 1 dataset : SRAFFAELE (paper assembly)
Vero - 4 datasets : KIM, DAVIDSON, TAIAROA, MATTHEWS

For each cell line we find the subsequent files:
*_plots_per_transcript.pdf -> this file contains the sharkfin plots for each transcript model of the assembly in which every point represent a position, its size is given by the number of datasets of that cell line in which the point is present (not significant, just present in the Nanocompore output). Points labeled with their letter sequence are those which have a LOR and a pvalue over the set threshold in all datasets of the cell line and that DO NOT fall in an IVT or ORF junction. The pvalue that they have in the plot is the max pvalue between the datasets and the LOR corresponding to it. Colors represent RAPID fragments.

*_plots_per_transcript_5p.pdf -> same as the *_plots_per_transcript.pdf but focused ONLY on sites with genomic Position <= 100. In this case the filter on junctions has been removed so labeled points MAY fall on junction sites.

*sites.txt -> this file contains the amount of modifications in a cell line per dataset and those which are shared. Columns are: ref_id (transcript model id), ORF (ORF encoded in that transcript), canonicity (if the transcript is canonical or not), shared_signif_mods (number of modification sites which are significant according to the established threshold in ALL the datasets of that cell line for that specific transcript model), shared_and_burrows ( number of modification sites which are significant according to the established threshold in ALL the datasets of that cell line for that specific transcript model and which are modified in Burrows paper), columns with sample names (number of modifications significant in that specific dataset for that transcript model).

*sites_identity.txt -> this file contains the identity of the SHARED modifications among all datasets of the specific cell line and all the information regarding them. (NO JUNCTION SITES)

*sites_identity_5p.txt -> this file contains the identity of the SHARED modifications among all datasets of the specific cell line AT THE 5P OF THE GENOME and all the information regarding them. (MAY BE JUNCTION SITES AND GENOMIC POS <=100)

*sites_identity_burrows.txt -> this file contains the identity of modifications sites listed in burrows paper with the relative max GMM pvalue and LOR. The column "shared" counts in how many datasets of that cell line, that site is considered significant for the chosen thresholds


The directory "shared" contains the analysis across all cell lines:

site_identity.txt -> this file is the equivalent of the files with the same name for each cell line.

site_identity_5p.txt -> this file is the equivalent of the files with the same name for each cell line. 

burrows_site_identity.txt -> this file is the equivalent of the files with the same name for each cell line.


PARAMETERS:
LOR_thresh <- 1
pval_thresh <- 0.01
IVT_junc_interval_left <- 25
IVT_junc_interval_right <- 25
ORF_junc_interval_left <- 15
ORF_junc_interval_right <- 15
burrows_sites <- c(22322,23317,27164,28417,28759,28927,29418)


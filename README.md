`[nicassiolab] (https://github.com/nicassiolab)`

# Nanopore analysis on SARS-CoV-2 pseudouridylation 

This directory contains scripts for the Nanopore analysis in "Discovering host protein interactions specific for SARS-CoV-2 RNA genome" by Giambruno et al. (doi: https://doi.org/10.1101/2022.07.18.499583). The analysis is performed on SARS-CoV-2 infected datasets from different cell line, with respect both to the viral SARS-CoV-2 genome or to the transcriptome NRCeq assembly (https://doi.org/10.1093/nar/gkac144). 

## How to run the analysis
Every submodule of the directory *scripts* contains a README.md file with the instructions to run the different parts of the analysis.
Edit the general configuration file present in the *general* directory before running the analysis.

## Content
* *files* : directory containing all the files necessary for the analysis
* *scripts* : directory containing all modules of the analysis and a general configuration file
  * *general* : directory containing a general configuration file to be set for the analysis
  * *alignments* : directory containing scripts to align reads to reference files
  * *gRNA_mods_detection* : directory containing scripts to detect modifications on genomic viral RNAs
  * *sgRNAs_mods_detection* : directory containing scripts to detect modifications on subgenomic viral RNAs


# gRNAs modifications detection

This directory contains scripts to detect putative modifications on genomic SARS-CoV-2 RNAs. 

## How to run the analysis

### Execute nanocompore analysis
To run the analysis compile the configuration file *configuration.sh* with the variables needed.
The analysis can be run just executing the script *gRNA_mods_detection.sh* as a bash script:
```bash
  ./gRNA_mods_detection.sh "$(pwd)" 
```
 or on a pbs cluster inserting the number of threads to use and the memory allocated for the job:
```bash
qsub -l select=1:ncpus=<threads>:mem=<memory to allocate> -v path="$(pwd)" gRNA_mods_detection.sh
```

### Generate results
To generate the results from the nanocompore analysis run the script *WT_vs_IVT_gRNAs_plots.R* on Rstudio. 

## Content
* *gRNA_mods_detection.sh* : script to run the modification analysis on genomic viral RNAs
* *config.sh* : configuration file with name of variables
* *images.sh* : script to download docker images and create image path
* *extract_reads_bam.py* : script to extract a list of reads from bamfile
* *WT_vs_IVT_gRNAs_plots.R* : script to generate results and WT_vs_IVT_gRNAs_plots
* *functions.R* : script containing functions used in *WT_vs_IVT_gRNAs_plots.R*
* *variables.R* : script containing variables used in *WT_vs_IVT_gRNAs_plots.R*

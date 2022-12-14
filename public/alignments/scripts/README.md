
# Alignments

This directory contains scripts to map SARS-CoV-2 infected datasets from different cell line to the viral SARS-CoV-2 genome or to the transcriptome NRCeq assembly (https://doi.org/10.1093/nar/gkac144). 

## How to run the analysis

### Execute nanocompore analysis
To run the analysis compile the configuration file *config.sh*  with the variables needed.
The analysis can be run just executing the script *per_cell_line_viral_alignments.sh* or *viral_alignments.sh* as a bash script or on a pbs cluster.
```bash
  ./per_cell_line_viral_alignments.sh
  ./viral_alignments.sh 
```

## Content
* *per_cell_line_viral_alignments.sh* : script that groups datasets based on the cell line and subsequently align them to the reference viral genome and transcriptome assembly.
* *viral_alignments.sh* : script that aligns each dataset to the reference viral genome and transcriptome assembly.
* *config.sh* : configuration file with name of variables
* *images.sh* : script to download docker images and create image path

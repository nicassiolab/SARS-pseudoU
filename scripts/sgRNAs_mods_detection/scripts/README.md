
# sgRNAs modifications detection

This directory contains scripts to detect putative modifications on sub-genomic SARS-CoV-2 RNAs using NRCeq transcriptome assembly from *Ugolini et al.* as a reference. 

## How to run the analysis

### Execute nanocompore analysis
To run the analysis compile the configuration file *configuration.sh* with the variables needed.
The analysis can be run just executing the script *nanocompore_per_cell_line_viral.sh* as a bash script:

```bash
  ./nanocompore_per_cell_line_viral.sh "$(pwd)" 
```
 or on a pbs cluster inserting the number of threads to use and the memory allocated for the job:
```bash
qsub -l select=1:ncpus=<threads>:mem=<memory to allocate> -v path="$(pwd)" nanocompore_per_cell_line_viral.sh 
```
This will generate Nanocompore results for every cell line analysed and will generate also violin plots for specific sites.


### Generate results
To generate the final results grouped per cell line from the nanocompore analysis run the script *sites_analysis/cell_line_per_cell_line_analysis.R* on Rstudio.
This analysis will produce an .xls file for each cell line whose results are described in the paper. 
These files will be the input for the next step of the analysis.
To generate final plots and results shared between different cell lines, run the script *sites_analysis/allfiles_all_cell_lines_analysis.R* on Rstudio. 


## Content
* *gRNA_mods_detection.sh* : script to run the modification analysis on genomic viral RNAs
* *config.sh* : configuration file with name of variables
* *images.sh* : script to download docker images and create image path
* *nanocompore_per_cell_line_viral.sh* : script to run the nanocompore analysis on data for each cell line
* *pseudoU_signal_plots.py* : python script to generate violin plots of electrical signal comparison for each k-mer
* *sites_analysis/functions.R* : script containing functions used in *sites_analysis/allfiles_all_cell_lines_analysis.R* and *sites_analysis/cell_line_per_cell_line_analysis.R*
* *sites_analysis/variables.R* : script containing variables used in *sites_analysis/allfiles_all_cell_lines_analysis.R* and *sites_analysis/cell_line_per_cell_line_analysis.R*
* *peakcalling/peakcalling.sh* : script to start pyhton peakcalling for every cell line analysed
* *peakcalling/nanocompore_peak_calling.py* : python script to perform peakcalling on close significant k-mer

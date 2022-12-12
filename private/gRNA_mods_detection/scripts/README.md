# Script directory for detecting modifications on genomic RNAs 

## PsiNanopore
Tool to identify mismatches between IVT control and WT sample. The script is edited with some modifications on the hypothesized minus strand processing. 
To run the script:
- copy the script into the Psi-Nanopore-main github directory;
- copy files into the data directory;
- edit the viral genome in order to have no backspace inside the sequence.
The execution gave 0 pseudoU sites.

## marginAlign
Tool to call variations in ONT reads. Tried running but failed in creation of virtualenv with suggested commands. Failed installation of pysam with python2.7 (tried also with docker image, conda and yum).

## nanoRMS
Tool to find pseudouridines in datasets. 
Results:
- viral alignments to genome: not working (error: OverflowError: Python int too large to convert to C long)


#!/bin/bash

# parameters to edit
BASECALLING="guppy_initial"				# indicate basecalling version (since the initial basecalling for the different datasets has been performed with different versions of guppy, we just call it "guppy_initial")
which_nucl="U"						# nucleotide for which identifying the modification
condition="WT"						# condition of processed samples (can be an array for multiple conditions)
THREADS=20						# 20 is indicated for f5c, while 3 is indicated for nanocompore 
PARALLEL_JOBS=5						# number of parallel jobs run by f5c eventalign

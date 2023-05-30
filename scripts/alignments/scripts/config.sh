#!/bin/bash

# parameters to edit
THREADS=20                                              # number of threads
BASECALLING="guppy_initial"                             # # indicate basecalling version (since the initial basecalling for the different datasets has been performed with different versions of guppy, we just call it "guppy_initial")
SAMPLE_CONDITION="WT"                                   # list of condition of the samples to be processed
condition_per_cell_line="WT"				# condition of samples to be processed collectively per cell line
PARALLEL_JOBS=5						# numner of parallel jobs for f5c eventalign

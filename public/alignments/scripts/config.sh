#!/bin/bash

# parameters to edit
THREADS=10                                              # number of threads
BASECALLING="guppy_initial"                             # version of guppy used for basecalling
SAMPLE_CONDITION="WT"                                   # list of condition of the samples to be processed
condition_per_cell_line="WT"				# condition of samples to be processed collectively per cell line
PARALLEL_JOBS=5						# numner of parallel jobs for f5c eventalign

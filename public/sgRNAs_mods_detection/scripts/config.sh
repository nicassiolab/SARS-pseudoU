#!/bin/bash

# parameters to edit
BASECALLING="guppy_initial"
condition="WT"						# condition of processed samples (can be an array for multiple conditions)
THREADS=20						# 20 is indicated for f5c, while 3 is indicated for nanocompore 
PARALLEL_JOBS=5						# number of parallel jobs run by f5c eventalign

library(dplyr)
library(tidyverse)
library(readr)

# script to compare T->C mismatches for IVT and gRNAs reads 

IVT_mismatches_file <- "/Volumes/scratch/TSSM/cugolini/cov/analysis/gRNA_mods/guppy_initial/pysamstats/IVT_T_C_mismatches.txt"
gRNAS_mismatches_file <-"/Volumes/scratch/TSSM/cugolini/cov/analysis/gRNA_mods/guppy_initial/pysamstats/grnas_T_C_mismatches.txt"

mismatches_IVT <-
  read_tsv(
    IVT_mismatches_file,
    col_names = c(
      "chrom",
      "genomicPos",
      "ref_nucl",
      "reads_all",
      "C_reads",
      "mism_percentage"
    ),
    col_types = c('cncnnn')
  )%>%
  mutate(genomicPos=(genomicPos-1))

mismatches_gRNAs <-
  read_tsv(
    gRNAS_mismatches_file,
    col_names = c(
      "chrom",
      "genomicPos",
      "ref_nucl",
      "reads_all",
      "C_reads",
      "mism_percentage"
    ),
    col_types = c('cncnnn')
  ) %>%
  mutate(genomicPos=(genomicPos-1))

# get overlap of pysam identified mismatches with gRNAs nanocompore k-mers
T_C_mism <- mismatches_gRNAs %>%
  filter(!genomicPos %in% mismatches_IVT$genomicPos)
mism_diff <- left_join(mismatches_gRNAs,mismatches_IVT, by=c("genomicPos","chrom","ref_nucl","reads_all","C_reads")) %>%
  subset(abs(mism_percentage.x-mism_percentage.y)>=0.01)

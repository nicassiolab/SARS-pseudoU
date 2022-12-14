#!/bin/bash
BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
WD="${BASEDIR}/analysis/sylamer"
NANOCOMPORE_FULL="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/caco2_sylamer_input.txt"
#sylamer=$BASEDIR/profiles/METTL3_polyA/bin/sylamer-12-342/sylamer

mkdir -p $WD

awk 'BEGIN{OFS=FS="\t"}NR>1 && $4<0.5{$4=-log($4)/log(10); print $1,$2,$3,$4}' $NANOCOMPORE_FULL | sort -k1,1 -k2,2 -k4,4gr | sort -k1,1 -k2,2n -u | sort -k4,4gr | awk '{print ">"$1"_"$2"\n"$3}' > full_sites_ext.fa
grep ">" full_sites_ext.fa | tr -d '>' > full_sites_ext.fa.ids.txt

#$sylamer -fasta $ANALYSIS/motifs/full_sites_ext.fa -universe $ANALYSIS/motifs/full_sites_ext.fa.ids.txt -k 5 -grow 100 --over  -o $ANALYSIS/motifs/sylamer_full.out
#$sylamer -fasta $ANALYSIS/motifs/full_sites_ext.fa -universe $ANALYSIS/motifs/full_sites_ext.fa.ids.txt -k 5 -grow 100 --logfold  -o $ANALYSIS/motifs/sylamer_full.lfc.out

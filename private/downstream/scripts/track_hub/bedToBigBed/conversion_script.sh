#!/bin/bash

#DATADIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/shared/bedtracks"
#WD="${DATADIR}/hub"
SCRIPTDIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/scripts/track_hub/bedToBigBed"

###     overlap over 3 databases

#wget https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.chrom.sizes
#mv $SCRIPTDIR/wuhCor1.chrom.sizes $WD
source activate $HOME/miniconda/envs/nanopore_pipeline_env/
####	 convert shared tracks
#mkdir -p $WD/bb
#for filename in $DATADIR/*.bed; do
#	name=${filename##*/}
#	base=${name%.bed}
#	sort -k1,1 -k2,2n $filename > $WD/"$base".sorted.bed      
#	$SCRIPTDIR/bedToBigBed $WD/"$base".sorted.bed $WD/wuhCor1.chrom.sizes $WD/bb/"$base".bigBed
#done


###	convert tracks for single cell lines
DATADIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/BEDTRACKS"
WD="${DATADIR}/hub"
mkdir -p $WD/bb

for filename in $DATADIR/*.bed; do
        name=${filename##*/}
        base=${name%.bed}
        sort -k1,1 -k2,2n $filename > $WD/"$base".sorted.bed
        $SCRIPTDIR/bedToBigBed $WD/"$base".sorted.bed /hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/results_allfiles_LOR05_pval001/shared/bedtracks/hub//wuhCor1.chrom.sizes $WD/bb/"$base".bigBed
done





#!/bin/bash
set -e -o pipefail
#export LC_ALL=C

source /hpcnfs/home/ieo5215/miniconda/etc/profile.d/conda.sh

SCRIPTDIR="/hpcnfs/scratch/FN/TL/cugolini/cov/scripts/downstream/scripts/peakcalling"
BASEDIR="/hpcnfs/scratch/TSSM/cugolini/cov"
DATA="$BASEDIR/data"
WD="$BASEDIR/analysis/per_cell_line"
REF_DATA="/hpcnfs/scratch/FN/camilla/nanopore/data"
ENVS="/hpcnfs/home/ieo5215/miniconda/envs"
nanocomp104="/hpcnfs/scratch/TSSM/cugolini/cov/img/nanocompore_v1.0.4.sif"

###	vero cells
mkdir -p  $WD/vero/NANOCOMPORE/peakcalling

for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/vero/NANOCOMPORE/sampcomp/*/outnanocompore_results.tsv; do
	base=${filename##*/sampcomp/}
	name=${base%/outnanocompore_results.tsv}
	singularity exec -B /hpcnfs/scratch/ $nanocomp104 python3 $SCRIPTDIR/nanocompore_peak_calling.py -i $filename -o $WD/vero/NANOCOMPORE/peakcalling/"$name".bed
done

###	caco2 cells
mkdir -p  $WD/caco2/NANOCOMPORE/peakcalling

for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/per_cell_line/caco2/NANOCOMPORE/sampcomp/*/outnanocompore_results.tsv; do
        base=${filename##*/sampcomp/}
        name=${base%/outnanocompore_results.tsv}
        singularity exec -B /hpcnfs/scratch/ $nanocomp104 python3 $SCRIPTDIR/nanocompore_peak_calling.py -i $filename -o $WD/caco2/NANOCOMPORE/peakcalling/"$name".bed
done


###	calu3 cells
mkdir -p  $WD/calu3/NANOCOMPORE/peakcalling

for filename in /hpcnfs/scratch/TSSM/cugolini/cov/analysis/7_samples_extraction/nanocompore/HUXELERATE_RESULTS/extraction/nanocompore/comparison/sraf_calu3/*/out_nanocompore_results.tsv; do
        base=${filename##*/sraf_calu3/}
        name=${base%/out_nanocompore_results.tsv}
        singularity exec -B /hpcnfs/scratch/ $nanocomp104 python3 $SCRIPTDIR/nanocompore_peak_calling.py -i $filename -o $WD/calu3/NANOCOMPORE/peakcalling/"$name".bed
done


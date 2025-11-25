#!/bin/bash

set -e -o pipefail

# load variables from general configuration file
if [ -z "$path" ]
then
      CURR_DIR=$1
else
      CURR_DIR=$path
fi

CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f7- |rev)                                    # obtain configuration file directory
source $CONFIG/scripts/scripts/general/config.sh


WD="$BASEDIR/analysis/sgRNAs_mods_detection/guppy_initial/structural_plots/TRSL"
template="$CURR_DIR/TRSL.sto"
tx="5820681b-6191-45dd-959c-96a8097aafdd|2255::NC_045512v2:14-29871"
nanocomp_res="$BASEDIR/analysis/sgRNAs_mods_detection/guppy_initial/CaCo2/nanocompore/sampcomp/5820681b-6191-45dd-959c-96a8097aafdd_2255_NC_045512v2_14-29871/outnanocompore_results.tsv"
real_start=0
rfam_id="TRSL"

# create conda environment
if [[ -d ${ENVS}/r2r ]]; then
        source activate ${ENVS}/r2r
else
        conda env create -f $FILES/r2r.yml -p ${ENVS}/r2r
	source activate ${ENVS}/r2r

r2r ${WD}/TRSL.meta ${WD}/TRSL.pdf
cat <(head -n -1 $template) <(python ${SCRIPTDIR2}/create_annotations.py $nanocomp_res $NANOCOMP_BED $template $tx 0.01 6 $real_start $rfam_id ${WD}/image_file.svg) <(echo //) > ${WD}/TRSL_annotated.sto
r2r ${WD}/TRSL.meta ${WD}/TRSL.pdf

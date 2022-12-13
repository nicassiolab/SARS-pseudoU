#!/bin/bash

set -eo pipefail

# load variables from general configuration file
CURR_DIR=$(dirname "$(realpath "$0")")                                                  # obtain current script directory
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)                                    # obtain configuration file directory
source $CONFIG/general/config.sh

# load local configuration file
source $CURR_DIR/config.sh
# load images
source $CURR_DIR/images.sh



# Obtain human transcriptome fasta

$SINGC bedparse gtf2bed $HG_GTF | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > $HG_BED  					 #from gtf(ensemble format to BED file)
$SINGC fastaFromBed -name -split -s -fi $HG_GENOME_FA -bed $HG_BED -fo - | perl -pe 's/>(.+)\(.\)$/>$1/' > $HG_DATA/transcriptome_fasta.fa      	 #get transcriptome fasta


#!/bin/bash
set -e -o pipefail

if [ -z "$path" ]
then
      CURR_DIR=$1
else
      CURR_DIR=$path
fi					                                                # obtain current script directory


CONFIG_GENERAL=$(echo $CURR_DIR | rev | cut -d'/' -f4- |rev)                            # obtain configuration file directory
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f2- |rev) 
source $CONFIG_GENERAL/general/config.sh

# load local configuration file
source $CONFIG/config.sh $CONFIG
# load images
source $CONFIG/images.sh $CONFIG


WD="$BASEDIR/analysis/analysis/sgRNAs_mods_detection/$BASECALLING"


for cell_line in CaLu3 CaCo2 VeroE6; do
	for filename in $WD/$cell_line/nanocompore/sampcomp; do
		base=${filename##*/sampcomp/}
		name=${base%/outnanocompore_results.tsv}
		mkdir -p $WD/$cell_line/nanocompore/peakcalling_"$which_nucl"
		$NANOCOMPORE python3 $CURR_DIR/nanocompore_peak_calling.py -i $filename -o $WD/$cell_line/nanocompore/peakcalling_"$which_nucl"/"$name".bed

	done
done


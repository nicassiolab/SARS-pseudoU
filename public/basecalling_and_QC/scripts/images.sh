#!/usr/bin/env bash

# load variables
CURR_DIR=$(dirname "$(realpath "$0")")                                                  # obtain current script directory
CONFIG=$(echo $CURR_DIR | rev | cut -d'/' -f3- |rev)                                    # obtain configuration file directory

source $CONFIG/general/config.sh                                                        # load variable configuration file

# pull image
if [ ! -f "$IMG/guppy_latest.sif" ]; then
        cd $IMG
        singularity pull docker://genomicpariscentre/guppy:latest
fi

# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/guppy_latest.sif"


# pull image
if [ ! -f "$IMG/pycoqc_2.5.2.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/pycoqc:2.5.2
fi
# singularity command
PYCOQC="singularity exec -B $MOUNT_DIR $IMG/pycoqc_2.5.2.sif"


#!/usr/bin/env bash

SCRIPTDIR="$(dirname "$(realpath "$0")")"

# load variables
source $SCRIPTDIR/config.sh


# pull image
if [ ! -f "$IMG/nrceq_pipeline_latest.sif" ]; then
        cd $IMG
        singularity pull docker://cugolini/nrceq_pipeline:latest
fi

# singularity command
SINGC="singularity exec -B $MOUNT_DIR $IMG/nrceq_pipeline_latest.sif"

# pull image for NanopolishComp
if [ ! -f "$IMG/nanocompore_latest.sif" ]; then
        cd $IMG
        singularity pull docker://adrienleger/nanocompore:latest
fi
# singularity command
NANOPOLISHCOMP="singularity exec -B $MOUNT_DIR $IMG/nanocompore_latest.sif"

# pull image for nanocompore
if [ ! -f "$IMG/nanocompore_v1.0.4.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/nanocompore:v1.0.4
fi
# singularity command
NANOCOMPORE="singularity exec -B $MOUNT_DIR $IMG/nanocompore_v1.0.4.sif"

# pull image for f5c
if [ ! -f "$IMG/f5c_v0.6.sif" ]; then
        cd $IMG
        singularity pull docker://tleonardi/f5c:v0.6
fi
# singularity command
F5C="singularity exec -B $MOUNT_DIR $IMG/f5c_v0.6.sif"

#!/bin/bash

# this script launches the singularity service for project combeffect

#FD_HOME=/gpfs/fs1/data/reddylab/Kuei/singularity/home
#FD_WORK=/gpfs/fs1/data/reddylab/Kuei
#FD_GITREPO=/gpfs/fs1/data/reddylab/Kuei/GitRepo
#FD_REDDYLAB=/gpfs/fs1/data/reddylab


SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${SCRIPT_DIR}/config_duke.sh

FD_HOME=${FD_SING}/home
N_PORT=9001

singularity exec \
    -H ${FD_HOME}:/home \
    -B ${FD_WORK}:/home/mount/work \
    -B ${FD_PRJ}:/home/mount/project \
    -B ${FD_RLAB}:/home/mount/reddylab \
    ${FD_SING}/singularity_proj_combeffect.sif \
    jupyter lab --NotebookApp.token="543@Psk" --no-browser --ip=0.0.0.0 --port=${N_PORT}


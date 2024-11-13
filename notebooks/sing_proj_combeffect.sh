#!/bin/bash

#source config_duke.sh

# this script launches the singularity service for project combeffect

#FD_HOME=/gpfs/fs1/data/reddylab/Kuei/singularity/home
#FD_WORK=/gpfs/fs1/data/reddylab/Kuei
#FD_GITREPO=/gpfs/fs1/data/reddylab/Kuei/GitRepo
#FD_REDDYLAB=/gpfs/fs1/data/reddylab

### Stackoverflow: how-can-i-get-the-source-directory-of-a-bash-script-from-within-the-script-itself
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${SCRIPT_DIR}/config_duke.sh

#singularity exec \
#    -H ${PWD}:/home \
#    -B ${FD_WORK}:/mount/work \
#    -B ${FD_PRJ}:/mount/project \
#    -B ${FD_RLAB}:/mount/reddylab \
#    ${FD_SING}/singularity_proj_combeffect.sif "$@"

#singularity exec \
#    -H ${PWD}:/home \
#    -B ${FD_WORK}:/home/mount/work \
#    -B ${FD_PRJ}:/home/mount/project \
#    -B ${FD_RLAB}:/home/mount/reddylab \
#    ${FD_SING}/singularity_proj_combeffect.sif "$@"
 
singularity exec \
    -H ${PWD}:/home \
    -B ${FD_WORK}:/mount/work \
    -B ${FD_PRJ}:/mount/project \
    -B ${FD_RLAB}:/mount/reddylab \
    ${FD_SING}/singularity_proj_combeffect.sif "$@"
    
#rm -rf ${PWD}/mount

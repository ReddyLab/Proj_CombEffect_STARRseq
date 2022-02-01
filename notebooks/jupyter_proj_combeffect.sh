#!/bin/bash
############################################################################
### this script launches the singularity service for project combeffect
############################################################################


### get flag ptions
### https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash
### https://stackoverflow.com/questions/16483119/an-example-of-how-to-use-getopts-in-bash
N_PORT=8888
usage() { echo "Usage: $0 [-p <PORT Number>]" 1>&2; exit 1; }
while getopts ':p:' flag; do
  case "${flag}" in
    p) 
        N_PORT=$OPTARG 
        echo "N_PORT: ${N_PORT}"
        ;;
    *) usage ;;
  esac
done
shift $((OPTIND -1))

### 
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source ${SCRIPT_DIR}/config_duke.sh
FD_HOME=${FD_SING}/home

### Run jupyter lab
singularity exec \
    -H ${FD_HOME}:/home \
    -B ${FD_WORK}:/mount/work \
    -B ${FD_PRJ}:/mount/project \
    -B ${FD_RLAB}:/mount/reddylab \
    ${FD_SING}/singularity_proj_combeffect.sif \
    jupyter lab --NotebookApp.token="543@Psk" --no-browser --ip=0.0.0.0 --port=${N_PORT}


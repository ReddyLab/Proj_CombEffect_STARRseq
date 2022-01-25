### set correct path based on the server I am at
FD_WORK=/home/mount/work
FD_RLAB=/home/mount/reddylab
FD_PRJ=/home/mount/project

### set working paths
FD_SRC=${FD_WORK}/source
FD_EXE=${FD_WORK}/exe
FD_ANN=${FD_WORK}/annotation

### set project related paths
FD_RES=${FD_WORK}/out/proj_combeffect
FD_LOG=${FD_RES}/log

### get flag ptions
### https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash
VERBOSE='false'
while getopts 'v' flag; do
  case "${flag}" in
    v) VERBOSE='true' ;;
    *) print_usage
       exit 1 ;;
  esac
done

### if verbose, print server and path
if ${VERBOSE}; then
    echo "You are in singularity_proj_combeffect"
    echo "BASE DIRECTORY:     ${FD_WORK}" 
    echo "PATH OF SOURCE:     ${FD_SRC}"
    echo "PATH OF EXECUTABLE: ${FD_EXE}"
    echo "PATH OF ANNOTATION: ${FD_ANN}"
    echo "PATH OF PROJECT:    ${FD_PRJ}"
    echo "PATH OF RESULTS:    ${FD_RES}"
    echo
fi

### load helper functions
source ${FD_PRJ}/config_func.sh
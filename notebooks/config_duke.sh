### Check which duke server I am at
### set correct path based on the server I am at
if echo $(pwd) | grep -q "gpfs"; then
    SERVER=HARDAC
    NODE=all
    
    FD_PREFIX="/gpfs/fs1"
    FD_WORK=${FD_PREFIX}/data/reddylab/Kuei
    FD_CODE=${FD_PREFIX}/data/reddylab/Kuei/GitRepo
    FD_RLAB=${FD_PREFIX}/data/reddylab
    FD_SING=${FD_WORK}/singularity
    
fi

if echo $(pwd) | grep -q "hpc"; then
    SERVER=DCC
    NODE=scavenger
    
    FD_PREFIX="/hpc"
    FD_WORK=/work/kk319
    FD_CODE=${FD_PREFIX}/home/kk319/GitRepo
    FD_RLAB=${FD_PREFIX}/group/reddylab
    FD_SING=${FD_RLAB}/Kuei/singularity
    
fi

### set working paths
FD_SRC=${FD_WORK}/source
FD_EXE=${FD_WORK}/exe
FD_ANN=${FD_WORK}/annotation

### set project related paths
FD_PRJ=${FD_CODE}/Proj_CombEffect_STARRseq/notebooks
FD_RES=${FD_WORK}/out/CombEffect_STARR
FD_LOG=${FD_RES}/log

### get flag ptions
### https://stackoverflow.com/questions/7069682/how-to-get-arguments-with-flags-in-bash
### https://stackoverflow.com/questions/16483119/an-example-of-how-to-use-getopts-in-bash
print_usage() { printf "Usage: hello"; }
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
    echo "You are on Duke Server: ${SERVER}"
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
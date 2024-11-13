### Check which duke server I am at
### set correct path based on the server I am at
if echo $(pwd -P) | grep -q "gpfs"; then
    SERVER=HARDAC
    NODE=all
    
    FD_PREFIX="/gpfs/fs1"
    FD_WORK=${FD_PREFIX}/data/reddylab/Kuei
    FD_CODE=${FD_PREFIX}/data/reddylab/Kuei/GitRepo
    FD_RLAB=${FD_PREFIX}/data/reddylab
    FD_SING=${FD_WORK}/singularity
    
    ### set working paths
    FD_ANN=${FD_WORK}/annotation
    FD_SRC=${FD_WORK}/source
    FD_EXE=${FD_WORK}/exe
    FD_BAC=${FD_WORK}/backup
fi

if echo $(pwd -P) | grep -q "hpc"; then
    SERVER=DCC
    NODE=scavenger
    
    FD_PREFIX="/hpc"
    FD_WORK=/work/kk319
    FD_CODE=${FD_PREFIX}/home/kk319/GitRepo
    FD_RLAB=${FD_PREFIX}/group/reddylab
    FD_SING=${FD_RLAB}/Kuei/singularity
    
    ### set working paths
    FD_ANN=${FD_RLAB}/Kuei/annotation
    FD_SRC=${FD_RLAB}/Kuei/source
    FD_EXE=${FD_RLAB}/Kuei/exe
fi

### set project related paths
FD_PRJ=${FD_CODE}/Proj_CombEffect_STARRseq/notebooks
FD_RES=${FD_WORK}/out/proj_combeffect
FD_LOG=${FD_RES}/log
FD_BIN=${FD_WORK}/bin

show_env() {
    echo "You are on Duke Server: ${SERVER}"
    echo "BASE DIRECTORY:     ${FD_WORK}" 
    echo "PATH OF SOURCE:     ${FD_SRC}"
    echo "PATH OF EXECUTABLE: ${FD_EXE}"
    echo "PATH OF ANNOTATION: ${FD_ANN}"
    echo "PATH OF PROJECT:    ${FD_PRJ}"
    echo "PATH OF RESULTS:    ${FD_RES}"
    echo
}
#fi

### load helper functions
source ${FD_PRJ}/config_func.sh

### container
FP_SIF=${FD_SING}/singularity_proj_combeffect.sif
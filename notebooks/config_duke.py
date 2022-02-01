import os
cwd = os.getcwd()

if "gpfs" in cwd:
    SERVER="HARDAC"
    
    FD_PREFIX="/gpfs/fs1"
    FD_WORK=os.path.join(FD_PREFIX, "data", "reddylab", "Kuei")
    FD_CODE=os.path.join(FD_PREFIX, "data", "reddylab", "Kuei", "GitRepo")
    FD_RLAB=os.path.join(FD_PREFIX, "data", "reddylab")
    FD_SING=os.path.join(FD_WORK,   "singularity")
    
    ### set working paths
    FD_ANN=os.path.join(FD_WORK, "annotation")
    FD_SRC=os.path.join(FD_WORK, "source")
    FD_EXE=os.path.join(FD_WORK, "exe")
    
if "hpc" in cwd:
    SERVER="DCC"
    
    FD_PREFIX="/hpc"
    FD_WORK="/work/kk319"
    FD_CODE=os.path.join(FD_PREFIX, "home/kk319/GitRepo")
    FD_RLAB=os.path.join(FD_PREFIX, "group/reddylab")
    FD_SING=os.path.join(FD_RLAB,   "Kuei/singularity")
    
    ### set working paths
    FD_ANN=os.path.join(FD_RLAB, "Kuei/annotation")
    FD_SRC=os.path.join(FD_RLAB, "Kuei/source")
    FD_EXE=os.path.join(FD_RLAB, "Kuei/exe")

### set project related paths
FD_PRJ=os.path.join(FD_CODE, "Proj_CombEffect_STARRseq/notebooks")
FD_RES=os.path.join(FD_WORK, "out/proj_combeffect")
FD_LOG=os.path.join(FD_RES,  "log")

if True:
    print(f"You are on Duke Server: {SERVER}")
    print(f"BASE DIRECTORY:     {FD_WORK}") 
    print(f"PATH OF SOURCE:     {FD_SRC}")
    print(f"PATH OF EXECUTABLE: {FD_EXE}")
    print(f"PATH OF ANNOTATION: {FD_ANN}")
    print(f"PATH OF PROJECT:    {FD_PRJ}")
    print(f"PATH OF RESULTS:    {FD_RES}")

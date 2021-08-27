### basic
library("tidyverse")
library("vroom")

### visualization
library("RColorBrewer")
library("gridExtra")

### set working path
get_base_path = function(){
    ### "/hpc/home/kk319"
    if (grepl("^/hpc", getwd())){
        cat("You are on Duke Server: DCC\n")
        return("/work/kk319")
    } 
    
    ### "/gpfs/fs1/data/reddylab/Kuei"
    if (grepl("/gpfs/fs1", getwd())){
        cat("You are on Duke Server: HARDAC\n")
        return("/data/reddylab/Kuei")
    }
    
    print("Server name not match DCC or HARDAC")
    return("Server name not match DCC or HARDAC")
}

FD_BASE = get_base_path()
FD_SRC  = file.path(FD_BASE, "source")
FD_EXE  = file.path(FD_BASE, "exe")
FD_ANN  = file.path(FD_BASE, "annotation")
FD_WRK  = file.path(FD_BASE, "out/CombEffect_STARR")
FD_LOG  = file.path(FD_WRK,  "log")

cat("BASE DIRECTORY:    ", FD_BASE, "\n")
cat("WORK DIRECTORY:    ", FD_WRK,  "\n")
cat("LOG  DIRECTORY:    ", FD_LOG,  "\n")
cat("PATH OF SOURCE:    ", FD_SRC,  "\n")
cat("PATH OF EXECUTABLE:", FD_EXE,  "\n")
cat("PATH OF ANNOTATION:", FD_ANN,  "\n")

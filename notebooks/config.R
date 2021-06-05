### basic
library("tidyverse")
library("vroom")
library("gtools")
library("corrr")

### visualization
library("pheatmap")
library("RColorBrewer")
library("gridExtra")
library("grid")

get_work_path = function(){
    ### "/hpc/home/kk319"
    if (grepl("^/hpc", getwd())){
        return("/work/kk319")
    } 
    
    ### "/gpfs/fs1/data/reddylab/Kuei"
    if (grepl("/gpfs/fs1", getwd())){
        return("/data/reddylab/Kuei")
    } 
    
    print("Server name not match DCC or HARDAC")
    return("Server name not match DCC or HARDAC")
}

FD_WRK = get_work_path()
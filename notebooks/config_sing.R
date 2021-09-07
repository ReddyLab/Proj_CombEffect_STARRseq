### basic
library("tidyverse")
library("vroom")

### visualization
library("RColorBrewer")
library("gridExtra")

### set paths
FD_WORK = "/home/mount/work"
FD_RLAB = "/home/mount/reddylab"
FD_PRJ  = "/home/mount/project"

FD_SRC  = file.path(FD_WORK, "source")
FD_EXE  = file.path(FD_WORK, "exe")
FD_ANN  = file.path(FD_WORK, "annotation")
FD_RES  = file.path(FD_WORK, "out", "CombEffect_STARR")
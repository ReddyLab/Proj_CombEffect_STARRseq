### basic
library("tidyverse")
library("vroom")
library("pryr")

### visualization
library("RColorBrewer")
library("cowplot")
library("gridExtra")
library("grid")

### set paths
FD_WORK = "/mount/work"
FD_RLAB = "/mount/reddylab"
FD_PRJ  = "/mount/project"

FD_SRC  = file.path(FD_WORK, "source")
FD_EXE  = file.path(FD_WORK, "exe")
FD_ANN  = file.path(FD_WORK, "annotation")
FD_RES  = file.path(FD_WORK, "out", "proj_combeffect")

cat("You are in singularity_proj_combeffect", "\n")
cat("BASE DIRECTORY:    ", FD_WORK, "\n")
cat("PATH OF SOURCE:    ", FD_SRC,  "\n")
cat("PATH OF EXECUTABLE:", FD_EXE,  "\n")
cat("PATH OF ANNOTATION:", FD_ANN,  "\n")
cat("PATH OF PROJECT:   ", FD_PRJ,  "\n")
cat("PATH OF RESULTS:   ", FD_RES,  "\n")
cat("")
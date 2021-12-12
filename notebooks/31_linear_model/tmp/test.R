### set env
source("/home/mount/project/config_sing.R")
cat("###################################################\n\n")

### Get argument: Chromomsome
ARGS        = commandArgs(trailingOnly=TRUE)
CHROM       = ARGS[1]
IS_INPUT20X = as.logical(ARGS[2])
FDIRY       = ARGS[3]

### set global variables
SAMPLES_TOT = c(
    paste0("Input", 1:5),
    paste0("Input", 1:5, "_20x"),
    paste0("TFX",   2:5, "_DMSO"),
    paste0("TFX",   2:5, "_Dex"))

SAMPLES_INP = c(
    paste0("Input", 1:5),
    paste0("TFX",   2:5, "_DMSO"),
    paste0("TFX",   2:5, "_Dex"))

SAMPLES_INP20X = c(
    paste0("Input", 1:5, "_20x"),
    paste0("TFX",   2:5, "_DMSO"),
    paste0("TFX",   2:5, "_Dex"))

if (IS_INPUT20X){
    SAMPLES = SAMPLES_INP20X
    FDIRY   = paste0(FDIRY, "_", "input20x")
} else {
    SAMPLES = SAMPLES_INP
}
FD_OUT = file.path(FD_RES, FDIRY)

### print start message
cat("Is Input20x used?", IS_INPUT20X, "\n")
cat("Output Directory:", FD_OUT, "\n")


### set environment
cat("\n++++++++++ Set environment  ++++++++++\n")
#source("/home/mount/project/config_sing.R")
source("config_sing.R")

cat("\n++++++++++ Set global variables ++++++++++\n")
### Get argument: Chromomsome
ARGS        = commandArgs(trailingOnly=TRUE)
CHROM       = ARGS[1]
IS_INPUT20X = as.logical(ARGS[2])
FDIRY       = ARGS[3]
THRESHOLD   = 10

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

if (IS_INPUT20X) {
    SAMPLES = SAMPLES_INP20X
    FDIRY   = paste0(FDIRY, "_", "input20x")
} else {
    SAMPLES = SAMPLES_INP
}
FD_OUT = file.path(FD_RES, "model_linear", FDIRY, CHROM)

### print start message
cat("Chromosome:      ", CHROM,       "\n")
cat("Is Input20x used?", IS_INPUT20X, "\n")
cat("Output Directory:", FD_OUT,      "\n")

###################################################
### Import library size
###################################################
cat("\n++++++++++ Import library size ++++++++++\n")

### Helper function to get
get_group = function(idn_sample){
    idn = idn_sample
    
    idn = str_replace(
        string = idn, 
        pattern = "Input[0-9]", 
        replacement = "Input")
    
    idn = str_remove(
        string = idn, 
        pattern = "_20x")
    
    idn = str_replace(
        string = idn, 
        pattern = "TFX[0-9]_", 
        replacement="TFX_")
    return(idn)
}

### set path
fdiry = file.path(FD_RES, "source")
fname = "library_size.txt"
fpath = file.path(fdiry, fname)

### import library size
ctypes = c(col_integer(), col_character())
cnames = c("Size", "FPath")
dat_lib = read_tsv(fpath, col_types=ctypes, col_names = cnames)

### remove the total size
dat_lib = dat_lib %>% dplyr::filter(FPath != "total")

### summarize info from the file path
### stackoverflow: Extract only folder name right before filename from full path
dat_lib = dat_lib %>% 
    #mutate(Sample = tools::file_path_sans_ext(basename(FPath))) %>%
    mutate(Sample = basename(dirname(FPath))) %>%
    mutate(Group = get_group(Sample))
dat_lib = dat_lib %>% dplyr::select(Size, Sample, Group)
print(dat_lib)

###################################################
### Import annotated fragments
###################################################
cat("\n++++++++++ Import annotated fragments ++++++++++\n")

### set column names and types
ctypes = c(col_character(), col_integer(), col_integer(), col_integer(),
           col_character(), col_integer(), col_integer(),
           col_character(), col_double(),  col_integer())
cnames = c("Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag",
           "Chrom_MTF",  "Start_MTF",  "End_MTF",
           "Motif", "Score", "Overlap")

### import bed files for each sample 
fdiry = file.path(FD_RES, "annotation_fragment", "filter_motif_score095")
fname = paste0(CHROM, ".bed.gz")

lst_dat = lapply(SAMPLES_TOT, function(sam){
    ### set path
    fpath = file.path(fdiry, sam, fname)
    print(fpath); flush.console()
    
    ### import data
    dat = vroom(fpath, col_types=ctypes, col_names=cnames) %>% mutate(Sample = sam)
    return(dat)
})

### arrange data
dat_ann_frag = bind_rows(lst_dat)
lst_dat      = NULL

### check environment
cat("Current memory used after import data:\n")
mem_used()


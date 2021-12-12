
###################################################
# Set environment
###################################################
cat("\n++++++++++ Set environment  ++++++++++\n")

#source("/home/mount/project/config_sing.R")
source("config_sing.R")

###################################################
# Set global variables
###################################################
cat("\n++++++++++ Set global variables ++++++++++\n")

### Get argument: Chromomsome
ARGS            = commandArgs(trailingOnly=TRUE)
TARGET          = as.character(ARGS[1])  # which chromosome or region to run
FDIRY           = as.character(ARGS[2])  # the name of the output folder
N_CORE          = as.integer(ARGS[3])    # number of cores to register during the parallelization
THRESHOLD_COVER = as.numeric(ARGS[4])    # threshold for the motif score filteration

### set global variables
SAMPLES_INP = c(
    paste0("Input", 1:5),
    paste0("TFX",   2:5, "_DMSO"),
    paste0("TFX",   2:5, "_Dex"))

SAMPLES_INP20X = c(
    paste0("Input", 1:5, "_20x"),
    paste0("TFX",   2:5, "_DMSO"),
    paste0("TFX",   2:5, "_Dex"))

LST_SAMPLES = list(input=SAMPLES_INP, input20x=SAMPLES_INP20X)

FD_OUT = file.path(FD_RES, "model_linear", FDIRY, TARGET)
dir.create(FD_OUT, recursive = TRUE, showWarnings = FALSE)

### set motifs
fdiry  = file.path(FD_RES, "model_linear", FDIRY, TARGET)
fname  = "count_*"
fglob  = file.path(fdiry, fname)
fpaths = Sys.glob(fglob)
MOTIFS = basename(fpaths)

### set column names and types
CTYPES = c(col_character(), col_character(), col_integer(), col_integer())
CNAMES = c("Sample", "Motif", "N_Motif", "Value")
           
### print start message
cat("Target:           ", TARGET,          "\n")
cat("Output Directory: ", FD_OUT,          "\n")
cat("#Cores Resgister: ", N_CORE,          "\n")
cat("Threshold (Cover):", THRESHOLD_COVER, "\n")

###################################################
# Import library size
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
dat_lib = read_tsv(fpath, col_types=ctypes, col_names=cnames)

### remove the total size
dat_lib = dat_lib %>% dplyr::filter(FPath != "total")

### summarize info from the file path
### stackoverflow: Extract only folder name right before filename from full path
dat_lib = dat_lib %>% 
    mutate(Sample = basename(dirname(FPath))) %>%
    mutate(Group = get_group(Sample))
dat_lib = dat_lib %>% dplyr::select(Size, Sample, Group)

###################################################
# Get motif count table
###################################################
cat("\n++++++++++ Get motif count table ++++++++++\n")

### start
registerDoParallel(N_CORE)
timer_start = Sys.time()

### loop through each motif to get the marginal effect
#lst_tmp = foreach(motif = MOTIFS) %do% {

#for (motif in MOTIFS) {
lst_tmp = foreach (motif = MOTIFS) %dopar% {
    ### start message and get the name of motif
    ### example: motif = "AHR_merge.bed.gz"
    mtf = str_remove_all(motif, pattern = "count_|\\.tsv")
    msg = paste(mtf, "Start")
    cat(msg, "\n"); flush.console()
    
    #++++++++++++++++++++++++++++++++++++++++++
    
    ### set directory
    fdiry  = file.path(FD_RES, "model_linear", FDIRY, TARGET)
    fname  = motif
    fpath  = file.path(fdiry, fname)
    
    ### PRINT
    msg = paste(mtf, "Import", fpath)
    cat(msg, "\n"); flush.console()
    
    ### import
    dat = read_tsv(fpath, show_col_types = FALSE)
    
    #++++++++++++++++++++++++++++++++++++++++++
    
    ### PRINT
    msg = paste(mtf, "Set")
    cat(msg, "\n"); flush.console()
    
    ### get total count fragments for each sample
    dat = dat %>% group_by(Sample) %>% summarize(Value = sum(Value))
    
    ### normalize counts by library size
    dat = dat %>% left_join(dat_lib, by="Sample")
    dat = dat %>%
        mutate(Norm_Value    = Value / Size) %>% 
        mutate(Lognorm_Value = log2(Value) - log2(Size)) %>%
        mutate(X = Group)
    
    #++++++++++++++++++++++++++++++++++++++++++
    
    lst_res = lapply(LST_SAMPLES, function(SAMPLES){
        ###
        tmp = dat %>% dplyr::filter(Sample %in% SAMPLES)

        ### create design matrix
        idxs  = c("Input", "TFX_DMSO", "TFX_Dex")
        tmp$X = factor(tmp$X, levels=idxs)
        X = model.matrix(~X, tmp)

        ### setup design matrix
        X[,"XTFX_DMSO"] = X[,"XTFX_DMSO"] + X[,"XTFX_Dex"]

        ### fit model and get the summary
        y = tmp$Norm_Value    
        fit = lm(y ~ X + 0)

        y = tmp$Lognorm_Value
        fit_log = lm(y ~ X + 0)

        ### arrange
        lst = list()
        lst$cnt     = tmp
        lst$fit     = fit
        lst$fit_log = fit_log
        lst$X       = X
        return(lst)
    })
    
    #++++++++++++++++++++++++++++++++++++++++++
    
    ### store the results
    fdiry = FD_OUT
    fname = paste0("lm_", mtf, ".RDS") # str_replace(mtf, pattern = "/", replacement = "_")
    fpath = file.path(fdiry, fname)
    
    ### PRINT
    msg = paste(mtf, "Store")
    cat(msg, "\n"); flush.console()
    
    ###
    saveRDS(lst_res, fpath)
    
} # end loop (foreach)

### PRINT
msg = "Done!"
cat(msg, "\n"); flush.console()

### print end message
timer = Sys.time()
cat("Timer of the loop:\n")
print(timer - timer_start)


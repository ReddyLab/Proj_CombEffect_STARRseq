
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
IS_INPUT20X     = as.logical(ARGS[2])    # is the new input being used
IS_LOG          = as.logical(ARGS[3])
FDIRY           = as.character(ARGS[4])  # the name of the output folder
N_CORE          = as.integer(ARGS[5])    # number of cores to register during the parallelization
THRESHOLD_COVER = as.integer(ARGS[6])    # threshold for the low coverage filteration
THRESHOLD_MOTIF = as.numeric(ARGS[7])    # threshold for the motif score filteration

### set global variables
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

if (IS_LOG) {
    FDIRY   = paste0(FDIRY, "_", "log")
} 

FD_OUT = file.path(FD_RES, "model_linear", FDIRY, TARGET)
dir.create(FD_OUT, recursive = TRUE, showWarnings = FALSE)

### set motifs
fdiry  = file.path(FD_RES, "annotation_fragment", SAMPLES[1], TARGET)
fname  = "*_merge.bed.gz"
fglob  = file.path(fdiry, fname)
fpaths = Sys.glob(fglob)
MOTIFS = basename(fpaths)

### set column names and types
CTYPES = c(col_character(), col_integer(), col_integer(), col_integer(),
           col_character(), col_integer(), col_integer(),
           col_character(), col_double(),  col_integer())
CNAMES = c("Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag",
           "Chrom_MTF",  "Start_MTF",  "End_MTF",
           "Motif", "Score", "Overlap")
           

### print start message
cat("Target:           ", TARGET,          "\n")
cat("Is Input20x used? ", IS_INPUT20X,     "\n")
cat("Is Log transform? ", IS_LOG,          "\n")
cat("Output Directory: ", FD_OUT,          "\n")
cat("#Cores Resgister: ", N_CORE,          "\n")
cat("Threshold (Cover):", THRESHOLD_COVER, "\n")
cat("Threshold (Motif):", THRESHOLD_MOTIF, "\n")

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
dat_lib = read_tsv(fpath, col_types=ctypes, col_names = cnames)

### remove the total size
dat_lib = dat_lib %>% dplyr::filter(FPath != "total")

### summarize info from the file path
### stackoverflow: Extract only folder name right before filename from full path
dat_lib = dat_lib %>% 
    mutate(Sample = basename(dirname(FPath))) %>%
    mutate(Group = get_group(Sample))
dat_lib = dat_lib %>% dplyr::select(Size, Sample, Group)

###################################################
# Run motif analysis
###################################################
cat("\n++++++++++ Run motif analysis ++++++++++\n")

### start
#registerDoParallel(N_CORE)
timer_start = Sys.time()

### loop through each motif to get the marginal effect
lst_res = foreach(fname = MOTIFS[1:10]) %do% {
    
    ### start message and get the name of motif
    mtf = str_remove_all(fname, pattern = "_merge.bed.gz")
    msg = paste(mtf, "Start")
    cat(msg, "\n"); flush.console()
    
    ### import fragment annotation
    fdiry  = file.path(FD_RES, "annotation_fragment")
    lst_dat = lapply(SAMPLES, function(sam){
        ### set path
        fpath = file.path(fdiry, sam, TARGET, fname)    
        msg = paste(mtf, "Import", fpath)
        cat(msg, "\n"); flush.console()
        
        ### import data
        dat = read_tsv(fpath, col_types=CTYPES, col_names=CNAMES)
        if (nrow(dat) == 0){
            return(NULL)
        } else {
            ###
            num1 = nrow(dat)    
            dat = dat %>% 
                mutate(Sample = sam) %>%
                mutate(Length_MTF = End_MTF - Start_MTF)  %>%
                mutate(Length_Dif = Length_MTF - Overlap) %>% 
                dplyr::filter(Length_Dif == 0) %>%
                dplyr::filter(Score >= THRESHOLD_MOTIF)
            num2 = nrow(dat)
            
            ###
            msg = paste(num1, num2, sep="-")
            msg = paste(mtf, "Filter", sam, msg)
            cat(msg, "\n"); flush.console()
        }
    })
    
    ### arrange data
    dat = bind_rows(lst_dat)
    msg = paste(mtf, "Total", nrow(dat))
    cat(msg, "\n"); flush.console()
    
    ###################################################
    # Preprocess
    ###################################################
    #cat("+++++ Preprocess +++++\n")
    
    ### Filter out empty data
    if(nrow(dat) == 0){
        msg = paste(mtf, "Skip Empty")
        cat(msg, "\n"); flush.console()
        return(msg)
    }
    
    ### Filter: fully cover the motif and motif score
    #num1 = nrow(dat)
    #dat = dat %>% 
    #    dplyr::filter(Length_Dif == 0) %>%
    #    dplyr::filter(Score >= THRESHOLD_MOTIF)
    #num2 = nrow(dat)
    #msg = paste(num1, num2, sep="-")
    #msg = paste(mtf, "Filter", msg)
    #cat(msg, "\n"); flush.console()
    
    ### Filter out empty data    
    #if(nrow(dat) == 0){
    #    msg = paste(mtf, "Filter Empty")
    #    cat(msg, "\n"); flush.console()
    #    return(msg)
    #}
    
    ### Filter: No/Low coverage
    cnt = sum(dat$Count_Frag)
    if(cnt <= THRESHOLD_COVER){
        msg = paste(mtf, "Filter Low_Coverage")
        cat(msg, "\n"); flush.console()
        return(msg)
    }
    
    ###################################################
    # Create Count Table
    ###################################################
    #cat("+++++ Create Count Table +++++\n")
    
    ### get fragments
    dat = dat %>% 
        dplyr::select(Chrom_Frag, Start_Frag, End_Frag, Count_Frag, Sample) %>%
        distinct()

    ### get count for each sample
    dat = dat %>% group_by(Sample) %>% summarise(Value = sum(Count_Frag))
    
    ### normalize counts by library size
    dat = dat %>% left_join(dat_lib, by="Sample")
    dat = dat %>%
        mutate(Norm_Value    = Value / Size) %>% 
        mutate(Lognorm_Value = log10(Value) - log10(Size)) %>%
        mutate(X = Group)
    
    ###################################################
    # Analyze w/ Linear Model
    ###################################################
    #cat("+++++ Analyze w/ Linear Model +++++\n")
    
    ### create design matrix
    idxs  = c("Input", "TFX_DMSO", "TFX_Dex")
    dat$X = factor(dat$X, levels=idxs)
    X = model.matrix(~X, dat)
    if (IS_LOG){
        y = dat$Lognorm_Value
    } else {
        y = dat$Norm_Value
    }

    ### setup design matrix
    X[,"XTFX_DMSO"] = X[,"XTFX_DMSO"] + X[,"XTFX_Dex"]
    
    ### fit model and get the summary
    fit = lm(y ~ X + 0)
    
    ### arrange
    lst = list()
    lst$fit = fit
    lst$cnt = dat
    lst$X   = X
    lst$y   = y
    
    ### store the results
    fdiry = FD_OUT
    fname = paste0(mtf, ".RDS") # str_replace(mtf, pattern = "/", replacement = "_")
    fpath = file.path(fdiry, fname)
    saveRDS(lst, fpath)
    
    ### end message
    msg = paste(mtf, "Done")
    cat(msg, "\n"); flush.console()
    return(msg)
}

### print end message
timer = Sys.time()
cat("Timer of the loop:\n")
print(timer - timer_start)


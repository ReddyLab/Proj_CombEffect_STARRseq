
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
cat("Target:           ", TARGET,    "\n")
cat("Is Input20x used? ", IS_INPUT20X,     "\n")
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
    #mutate(Sample = tools::file_path_sans_ext(basename(FPath))) %>%
    mutate(Sample = basename(dirname(FPath))) %>%
    mutate(Group = get_group(Sample))
dat_lib = dat_lib %>% dplyr::select(Size, Sample, Group)

###################################################
# Run motif analysis
###################################################
cat("\n++++++++++ Run motif analysis ++++++++++\n")

### combination of motifs
dat_comb = t(combn(MOTIFS[1:7], 2))
lst_motif_pair = split(dat_comb, seq(nrow(dat_comb)))

### start
registerDoParallel(N_CORE)
timer_start = Sys.time()

### loop through each pair of motif to get the interactive effect
lst_res = foreach(motif_pair = lst_motif_pair) %dopar% {
    
    ### start message
    mtfs = sapply(motif_pair, function(fname){
        mtf = str_remove_all(fname, pattern = "_merge.bed.gz")
        return(mtf)
    })
    
    msg_mtf = paste(mtfs, collapse=" ")
    msg_mtf = paste(msg_mtf, "|")
    msg     = paste(msg_mtf, "Start")
    cat(msg, "\n"); flush.console()
    
    ### import annotated fragments for each motif
    lst_dat = lapply(motif_pair, function(fname){
        lst = lapply(SAMPLES, function(sam){

            ### set path
            fdiry = file.path(FD_RES, "annotation_fragment")
            fpath = file.path(fdiry, sam, TARGET, fname)    

            ### import data
            dat = read_tsv(fpath, col_types=CTYPES, col_names=CNAMES)
            if (nrow(dat) == 0){
                return(NULL) # exception handling: empty file
            } else {
                dat = dat %>% 
                    mutate(Sample = sam) %>%
                    mutate(Length_MTF = End_MTF - Start_MTF) %>%
                    mutate(Length_Dif = Length_MTF - Overlap)
                return(dat)
            }
        })

        ### arrange data
        dat = bind_rows(lst)
        return(dat)
    })
    
    ###################################################
    # Preprocess
    ###################################################
    #cat("+++++ Preprocess +++++\n")
    
    ### apply the same preprocessing for each motif within a motif pair
    for (idx in seq_along(lst_dat)) {
        ### extract info
        dat = lst_dat[[idx]]
        mtf = mtfs[idx]
        
        ### Filter out empty data
        if(nrow(dat) == 0){
            msg = paste(msg_mtf, mtf, "Skip Empty")
            cat(msg, "\n"); flush.console()
            return(msg)
        }

        ### Filter: fully cover the motif and motif score
        num1 = nrow(dat)
        dat = dat %>% 
            dplyr::filter(Length_Dif == 0) %>%
            dplyr::filter(Score >= THRESHOLD_MOTIF)
        num2 = nrow(dat)
        msg = paste(num1, num2, sep="-")
        msg = paste(msg_mtf, mtf, "Filter", msg)
        cat(msg, "\n"); flush.console()

        ### Filter out empty data    
        if(nrow(dat) == 0){
            msg = paste(mtf, "Filter Empty")
            cat(msg, "\n"); flush.console()
            return(msg)
        }

        ### Filter: No/Low coverage
        cnt = sum(dat$Count_Frag)
        if(cnt <= THRESHOLD_COVER){
            msg = paste(msg_mtf, mtf, "Filter Low_Coverage")
            cat(msg, "\n"); flush.console()
            return(msg)
        }
    } # end inner for loop
    
    ### arrange data after preprocessing
    df1  = lst_dat[[1]]
    df2  = lst_dat[[2]]
    mtf1 = unique(df1$Motif)
    mtf2 = unique(df2$Motif)
    lst_dat = NULL
    
    ###################################################
    # Create Count Table
    ###################################################
    #cat("+++++ Create Count Table +++++\n")
    
    ### extract fragments
    dat1 = bind_rows(df1, df2) %>% 
        dplyr::select(Chrom_Frag, Start_Frag, End_Frag, Count_Frag, Sample) %>%
        distinct
    dat2 = df1 %>% 
        dplyr::select(Chrom_Frag, Start_Frag, End_Frag, Count_Frag, Sample, Motif)
    dat3 = df2 %>% 
        dplyr::select(Chrom_Frag, Start_Frag, End_Frag, Count_Frag, Sample, Motif)
    
    ### match fragments for the motif pair
    dat = dat1 %>%
        full_join(dat2, by = c("Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag", "Sample")) %>%
        full_join(dat3, by = c("Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag", "Sample")) %>%
        mutate(Motif = paste(Motif.x, Motif.y, sep = "_")) %>%
        mutate(Motif = str_remove(string=Motif, pattern="_NA|NA_")) %>% 
        mutate(Group = str_remove(string = Sample, pattern = "[0-9]")) %>%
        mutate(X     = paste(Group, Motif, sep="_")) %>%
        mutate(X     = ifelse(str_detect(X, "Input"), "Input", X))
    
    ### get count for each sample
    dat = dat %>% 
        group_by(Sample, X) %>% 
        summarise(Value = sum(Count_Frag), .groups = 'drop')
    
    ### normalize counts by library size
    dat = dat %>% left_join(dat_lib, by="Sample")
    dat = dat %>%
        mutate(Norm_Value    = Value / Size) %>% 
        mutate(Lognorm_Value = log10(Value) - log10(Size))
    
    ###################################################
    # Analyze w/ Linear Model
    ###################################################
    #cat("+++++ Analyze w/ Linear Model +++++\n")
    
    ### create design matrix
    idx11 = paste("TFX_DMSO", mtf1,       sep="_")
    idx12 = paste("TFX_DMSO", mtf2,       sep="_")
    idx13 = paste("TFX_DMSO", mtf1, mtf2, sep="_")
    idx21 = paste("TFX_Dex",  mtf1,       sep="_")
    idx22 = paste("TFX_Dex",  mtf2,       sep="_")
    idx23 = paste("TFX_Dex",  mtf1, mtf2, sep="_")
    idxs  = c("Input", idx11, idx12, idx21, idx22, idx13, idx23)
    dat$X = factor(dat$X, levels=idxs)
    X = model.matrix(~X, dat)
    if (IS_LOG){
        y = dat$Lognorm_Value
    } else {
        y = dat$Norm_Value
    }
    
    ### setup design matrix
    idx11 = paste("XTFX_DMSO", mtf1,       sep="_")
    idx12 = paste("XTFX_DMSO", mtf2,       sep="_")
    idx13 = paste("XTFX_DMSO", mtf1, mtf2, sep="_")
    idx21 = paste("XTFX_Dex",  mtf1,       sep="_")
    idx22 = paste("XTFX_Dex",  mtf2,       sep="_")
    idx23 = paste("XTFX_Dex",  mtf1, mtf2, sep="_")
    X[,idx11] = X[,idx11] + X[,idx13] + X[,idx21] + X[,idx23]
    X[,idx12] = X[,idx12] + X[,idx13] + X[,idx22] + X[,idx23]
    X[,idx21] = X[,idx21] + X[,idx23]
    X[,idx22] = X[,idx22] + X[,idx23]
    X[,idx13] = X[,idx13] + X[,idx23]
    
    ### fit model and reduce the memory size
    fit = lm(y ~ X + 0)
    
    ### arrange the results
    lst = list()
    lst$fit = fit
    lst$cnt = dat
    lst$X   = X
    lst$y   = y
    
    ### store the results
    mtf1 = str_remove(string=motif_pair[1], pattern="_merge.bed.gz")
    mtf2 = str_remove(string=motif_pair[2], pattern="_merge.bed.gz")
    
    fdiry = FD_OUT
    fname = paste0(mtf1, "_", mtf2, ".RDS")
    fpath = file.path(fdiry, fname)
    dir.create(fdiry, recursive = TRUE, showWarnings = FALSE)
    saveRDS(lst, fpath)
    
    ### end message
    msg = paste(msg_mtf, "Done")
    cat(msg, "\n"); flush.console()
    return(msg)
}

### print end message
timer = Sys.time()
cat("Timer of the loop:\n")
print(timer - timer_start)


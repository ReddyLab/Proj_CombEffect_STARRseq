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
N_CORE      = as.integer(ARGS[4])
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
cat("#Cores Resgister:", N_CORE,      "\n")

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
    dat = read_tsv(fpath, col_types=ctypes, col_names=cnames) %>% mutate(Sample = sam)
    return(dat)
})

### arrange data
dat_ann_frag = bind_rows(lst_dat)
lst_dat      = NULL

### check environment
cat("Current memory used after import data:\n")
mem_used()

###################################################
### Preprocess
###################################################
cat("\n++++++++++ Preprocess ++++++++++\n")

### Filter: fully cover the motif
dat = dat_ann_frag
dat = dat %>% 
    mutate(Length_MTF = End_MTF - Start_MTF) %>%
    mutate(Length_Dif = Length_MTF - Overlap)
dat_ann_frag = NULL # release memory

cat("Filtering: fully cover the motif\n")
cat("    Before Filter:", "#Motif =", length(unique(dat$Motif)), "#Annot =", nrow(dat), "\n")
dat = dat %>% dplyr::filter(Length_Dif == 0)
cat("    After  Filter:", "#Motif =", length(unique(dat$Motif)), "#Annot =", nrow(dat), "\n")

### grouped by motif cluster and split the annotated fragments into list
lst    = dat %>% group_by(Motif) %>% group_split
motifs = lapply(lst, function(x){unique(x$Motif)}) %>% unlist
names(lst) = motifs
dat = NULL # release memory

lst_frag = lapply(lst, function(dat){
    tmp = dat %>% 
        group_by(Chrom_Frag, Start_Frag, End_Frag, Count_Frag, Motif, Sample) %>%
        summarize(N_Motif = n(), .groups = 'drop')
    return(tmp)
})

### filter out motifs that have almost no/low fragment in total
cat("Filtering: filter out motifs with no/low fragments\n")
cat("    Threshold:", THRESHOLD, "\n")
cat("    Before Filter:", "#Motif =", length(lst_frag), "\n")

lst = lst_frag
cnt = lapply(lst, function(dat){sum(dat$Count_Frag)})
lst = lst[cnt > 10]
lst_frag = lst

cat("    After  Filter:", "#Motif =", length(lst_frag), "\n")

### check environment
cat("Current memory used after import data:\n")
mem_used()

###################################################
### Set up linear model
###################################################
cat("\n++++++++++ Linear Model ++++++++++\n")

### Helper function
### https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/
stripGlmLR = function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()

  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  return(cm)
}

### start timer
timer_start = Sys.time()

### run linear model
registerDoParallel(N_CORE)
lst_tmp = lst_frag #head(lst_frag, 10)
lst_tmp = foreach(idn = names(lst_tmp)) %dopar% {

    ### extract and get fragments
    dat = lst_tmp[[idn]]
    dat = dat %>% 
        dplyr::select(Chrom_Frag, Start_Frag, End_Frag, Count_Frag, Sample) %>%
        distinct()

    ### count fragments for each sample
    dat = dat %>% group_by(Sample) %>% summarise(Value = sum(Count_Frag))
    
    ### choose samples for modeling
    dat = dat %>% dplyr::filter(Sample %in% SAMPLES)
    
    ### normalize counts by library size
    dat = dat %>% left_join(dat_lib, by="Sample")
    dat = dat %>%
        mutate(Norm_Value    = Value / Size) %>% 
        mutate(Lognorm_Value = log10(Value) - log10(Size)) %>%
        mutate(X = Group)

    ### create design matrix
    idxs  = c("Input", "TFX_DMSO", "TFX_Dex")
    dat$X = factor(dat$X, levels=idxs)
    X = model.matrix(~X, dat)
    y = dat$Norm_Value

    ### setup design matrix
    X[,"XTFX_DMSO"] = X[,"XTFX_DMSO"] + X[,"XTFX_Dex"]
    
    ### fit model and reduce the memory size
    fit = lm(y ~ X + 0)
    #fit = stripGlmLR(fit)
    res = summary(fit)
    res = stripGlmLR(res)
    
    ### arrange results
    lst = list()
    lst$res = res
    lst$cnt = dat
    lst$X   = X
    lst$y   = y
    
    ### store the results
    mtf   = str_replace_all(idn, pattern = "/", replacement = "_")
    #fdiry = file.path(FD_RES, "model_linear", FD_OUT, CHROM)
    fdiry = FD_OUT
    fname = paste0("motif_", mtf, ".RDS")
    fpath = file.path(fdiry, fname)
    
    #print(c(idn, mtf)); flush.console()
    dir.create(fdiry, recursive = TRUE, showWarnings = FALSE)
    saveRDS(lst, fpath)
}

### print end message
timer = Sys.time()
cat("Timer of the loop:\n")
print(timer - timer_start)


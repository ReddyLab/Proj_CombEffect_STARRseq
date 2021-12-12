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
THRESHOLD_COVER = as.integer(ARGS[2])    # threshold for the low coverage filteration
THRESHOLD_MOTIF = as.numeric(ARGS[3])    # threshold for the motif score filteration

SAMPLES_INP20X = c(
    paste0("Input", 1:5, "_20x"),
    paste0("TFX",   2:5, "_DMSO"),
    paste0("TFX",   2:5, "_Dex"))
SAMPLES = SAMPLES_INP20X

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
cat("Threshold (Cover):", THRESHOLD_COVER, "\n")
cat("Threshold (Motif):", THRESHOLD_MOTIF, "\n")

###################################################
# Set global variables
###################################################

### start
#registerDoParallel(N_CORE)
timer_start = Sys.time()

### loop through each motif to get the marginal effect
lst_res = foreach(fname = MOTIFS[1:5]) %do% {
    
    ### start loop timer
    timer = Sys.time()
    
    ### start message and get the name of motif
    mtf = str_remove_all(fname, pattern = "_merge.bed.gz")
    msg = paste(mtf, "Start")
    cat(msg, "\n"); flush.console()
    
    ### import fragment annotation
    fdiry  = file.path(FD_RES, "annotation_fragment")
    lst_dat = lapply(SAMPLES, function(sam){
        ### set path
        fpath = file.path(fdiry, sam, TARGET, fname)    
        
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
            return(dat)
        }
    })
    
    ### arrange data
    dat = bind_rows(lst_dat)
    
    ### end message
    msg = paste(mtf, "Done;", "nrow:", nrow(dat))
    cat(msg, "\n"); flush.console()
    print(Sys.time() - timer)
    return(msg)
}

### print end message
timer = Sys.time()
cat("Timer of the loop:\n")
print(timer - timer_start)



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
THRESHOLD_MOTIF = as.numeric(ARGS[4])    # threshold for the motif score filteration

### set global variables
SAMPLES = c(
    paste0("Input", 1:5),
    paste0("Input", 1:5, "_20x"),
    paste0("TFX",   2:5, "_DMSO"),
    paste0("TFX",   2:5, "_Dex"))

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
cat("Output Directory: ", FD_OUT,          "\n")
cat("#Cores Resgister: ", N_CORE,          "\n")
cat("Threshold (Motif):", THRESHOLD_MOTIF, "\n")

###################################################
# Set motif list
###################################################

### helper function
fun_chunk = function(x, n){ 
    if (n==1){ 
        ### EXCEPTION: split to only one chunk
        lst = list(x) 
    } else {
        ### split a vector into several chunks
        lst = split(x, cut(seq_along(x), n, labels = FALSE))
    }
    return(lst)
}

### split motifs into several chunks for parallel programming
#lst_motifs = fun_chunk(MOTIFS, N_CORE)

### combination of motifs
dat_comb = t(combn(MOTIFS, 2))
lst_motif_pairs = split(dat_comb, seq(nrow(dat_comb)))
#lst_motif_chunk = fun_chunk(lst_motif_pairs, N_CORE)

###################################################
# Get motif count table
###################################################
cat("\n++++++++++ Get motif count table ++++++++++\n")

### PRINT: start message
timer_start = Sys.time()

#registerDoParallel(cores=N_CORE)
#cl <- parallel::makeCluster(N_CORE)
#doParallel::registerDoParallel(cl)
#cl <- makeForkCluster(N_CORE)
#registerDoParallel(cl)

### loop through each motif to get the marginal effect
#foreach(index = 1:N_CORE) %do% {
    
    ### init
    #lst_motif_pair = lst_motif_chunk[[idx]]

### loop through each motif pair to estimate motif interaction
for (motif_pair in lst_motif_pairs){

    ### init
    is_created = FALSE
    mtfs = sapply(motif_pair, function(motif){
        mtf = str_remove_all(motif, pattern = "_merge.bed.gz")
        return(mtf)
    })
    mtfs = paste(mtfs, collapse="_")
    msg  = paste(mtfs, "Start")
    cat(msg, "\n"); flush.console()

    ###
    for (sam in SAMPLES){

         ### import annotated fragments for each motif
         lst_dat = lapply(motif_pair, function(motif){

             ###################################################
             # Import fragment annotation
             ###################################################

             ### SET: file path of annotated fragment
             fdiry = file.path(FD_RES, "annotation_fragment", sam, TARGET)
             fname = motif
             fpath = file.path(fdiry, fname)    

             ### PRINT: ready to import
             msg = paste(mtfs, sam, "Import", fpath)
             cat(msg, "\n"); flush.console()

             ### import data
             dat = read_tsv(fpath, col_types=CTYPES, col_names=CNAMES)

             ### HANDLE EXCEPTION: empty data
             if (nrow(dat) == 0){
                 msg = paste(mtfs, sam, "Skip Import_Empty")
                 cat(msg, "\n"); flush.console()
                 return(NULL)
             }

             ###################################################
             # Preprocess
             ###################################################

             ### FILTER:
             ###     filter out annotation not fully cover motif
             ###     filter out motif score lower than threshold
             num1 = nrow(dat)    
             dat = dat %>% 
                 mutate(Sample = sam) %>%
                 mutate(Length_MTF = End_MTF - Start_MTF)  %>%
                 mutate(Length_Dif = Length_MTF - Overlap) %>% 
                 dplyr::filter(Length_Dif == 0) %>%
                 dplyr::filter(Score >= THRESHOLD_MOTIF)
             num2 = nrow(dat)

             ### PRINT: result of filtering
             msg = paste(num1, num2, sep="-")
             msg = paste(mtfs, sam, "Filter", msg)
             cat(msg, "\n"); flush.console()

             ### HANDLE EXCEPTION: empty data after filteration
             if(nrow(dat) == 0){
                 msg = paste(mtfs, sam, "Skip Filter_Empty")
                 cat(msg, "\n"); flush.console()
                 return(NULL)
             }
             
             return(dat)
        }) # end lapply

        ### arrange data after preprocessing
        df1  = lst_dat[[1]]
        df2  = lst_dat[[2]]
        lst_dat = NULL

        ### HANDLE EXCEPTION: skip if one data is empty
        if (is.null(df1)){next}
        if (is.null(df2)){next}
        

        ###################################################
        # Create Count Table
        ###################################################
        cat("+++++ Create Count Table +++++\n")
        
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
            group_by(Sample, Motif.x, Motif.y, X) %>% 
            summarise(Value = sum(Count_Frag), .groups = 'drop')
        
        ###################################################
        # Store output
        ###################################################

        ### SET: file path for output count table
        fdiry = FD_OUT 
        fname = paste0("count_",  mtfs, ".tsv")
        fpath = file.path(fdiry, fname)

        ### store results
        ### create the table for the first sample or if the file is not yet created
        ### Otherwise, append the counts in the file
        if (is_created) {

            ### PRINT: file path for output count
            msg = paste(mtfs, sam, "Store_Append", fpath)
            cat(msg, "\n"); flush.console()

            ### append the file
            write.table(
                dat,
                file      = fpath,
                append    = TRUE,
                quote     = FALSE,
                sep       = "\t",
                row.names = FALSE,
                col.names = FALSE)

        } else {

            ### PRINT: file path for output count
            msg = paste(mtfs, sam, "Store_Create", fpath)
            cat(msg, "\n"); flush.console()

            ### create the file
            write.table(
                dat,
                file      = fpath,
                quote     = FALSE,
                sep       = "\t",
                row.names = FALSE,
                col.names = TRUE)

            ### update flag
            is_created = TRUE

        } # end if-else
        
    } # end inner for loop (SAMPLES)
} # end outer for loop (lst_motif_pair)
    
    

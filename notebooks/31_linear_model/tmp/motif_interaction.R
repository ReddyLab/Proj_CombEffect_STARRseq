### set env
source("/home/mount/project/config_sing.R")
cat("###################################################\n\n")

### Get argument: Chromomsome
ARGS  = commandArgs(trailingOnly=TRUE)
CHROM = ARGS[1]
print(CHROM)

###################################################
# Import library size
###################################################
cat("+++++ Import library size +++++\n")

### Helper function to get
get_sample = function(idn_sample){
    idn = idn_sample
    idn = str_replace(
        string = idn, 
        pattern = "Input[0-9]", 
        replacement = "Input")
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
cnames = c("Size", "Fpath")
dat_lib = read_tsv(fpath, col_types=ctypes, col_names = cnames)
dat_lib = dat_lib %>% 
    mutate(Sample = tools::file_path_sans_ext(basename(Fpath))) %>%
    mutate(Group = get_sample(Sample))
dat_lib = dat_lib %>% dplyr::select(Size, Sample, Group)

###################################################
# Import annotated fragments
###################################################
cat("+++++ Import annotated fragments +++++\n")

### set column names and types
ctypes = c(col_character(), col_integer(), col_integer(), col_integer(),
           col_character(), col_integer(), col_integer(),
           col_character(), col_double(),  col_integer())
cnames = c("Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag",
           "Chrom_MTF",  "Start_MTF",  "End_MTF",
           "Motif", "Score", "Overlap")

### set samples
SAMPLES = c(
    paste0("Input", 1:5),
    paste0("TFX",   2:5, "_DMSO"),
    paste0("TFX",   2:5, "_Dex"))

### import bed files for each sample 
fdiry = file.path(FD_RES, "annotation_fragment", "filter_motif_score095")
fname = paste0(CHROM, ".bed", ".gz") 

lst_dat = lapply(SAMPLES, function(sam){
    ### set path
    fpath = file.path(fdiry, sam, fname)
    print(fpath)
    
    ### import data
    dat = read_tsv(fpath, col_types=ctypes, col_names=cnames) %>% mutate(Sample = sam)
    return(dat)
})

### arrange data & clean memory
dat_ann_frag = bind_rows(lst_dat)
lst_dat = NULL

### check point: print message
head(dat_ann_frag, 3)
print(colnames(dat_ann_frag))

print("Dimension")
print(dim(dat_ann_frag))

print("Size")
print(object.size(dat_ann_frag), units="Mb")

###################################################
# Preprocess
###################################################

### grouped by motif cluster and 
### split the annotated fragments into list
cat("+++++ Preprocess: get motifs +++++\n")
dat    = dat_ann_frag
lst    = dat %>% group_by(Motif) %>% group_split
motifs = lapply(lst, function(dat){unique(dat$Motif)}) %>% unlist
names(lst) = motifs

### get the list
cat("+++++ Preprocess: get lst_frag +++++\n")
lst_frag = lapply(lst, function(dat){
    tmp = dat %>% 
        dplyr::select(Chrom_Frag, Start_Frag, End_Frag, Count_Frag, Motif, Sample) %>%
        distinct()
    return(tmp)
})

### filter out motifs that have almost no fragment in total
cat("+++++ Preprocess: filter motifs +++++\n")
cat("Before filteration: #Motifs =", length(lst_frag), "\n")

THRESHOLD = 10
lst = lst_frag
cnt = lapply(lst, function(dat){sum(dat$Count_Frag)})
lst = lst[cnt > 10]
lst_frag = lst

cat("Threshold =", THRESHOLD, "\n")
cat("After filteration: #Motifs =", length(lst_frag), "\n")

###################################################
# Set all pairs of motifs
###################################################
cat("+++++ Set all pairs of motifs +++++\n")

### combination of motifs
dat_comb = t(combn(names(lst_frag), 2))

### convert motif pairs into a list
lst_motif_pair = split(dat_comb, seq(nrow(dat_comb)))
names(lst_motif_pair) = lapply(
    lst_motif_pair, function(x){
        mtf1 = x[1]
        mtf2 = x[2]
        return(paste(mtf1, mtf2, sep="|"))
    } # end fun
) # end lapply

###################################################
# Get interactive effects
###################################################
cat("+++++ Get interactive effects +++++\n")

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

### loop through each pair of motifs
### estimate interaction effect of each motif pair
lst_res = lapply(lst_motif_pair, function(x){

    ### extract fragments for each motif
    mtf1 = x[1]
    mtf2 = x[2]
    df1 = lst_frag[[mtf1]]
    df2 = lst_frag[[mtf2]]
    
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
    
    ### annotate fragments based on motif annotation
    idx11 = paste("TFX_DMSO", mtf1,       sep="_")
    idx12 = paste("TFX_DMSO", mtf2,       sep="_")
    idx13 = paste("TFX_DMSO", mtf1, mtf2, sep="_")
    idx21 = paste("TFX_Dex",  mtf1,       sep="_")
    idx22 = paste("TFX_Dex",  mtf2,       sep="_")
    idx23 = paste("TFX_Dex",  mtf1, mtf2, sep="_")
    idxs  = c("Input", idx11, idx12, idx21, idx22, idx13, idx23)
    tmp = dat %>% 
        group_by(Sample, X) %>% 
        summarise(Value = sum(Count_Frag), .groups = 'drop')
    
    ### normalize counts by library size
    tmp = tmp %>% left_join(dat_lib, by="Sample")
    tmp = tmp %>%
        mutate(Norm_Value    = Value / Size) %>% 
        mutate(Lognorm_Value = log10(Value) - log10(Size))
    tmp$X = factor(tmp$X, levels=idxs)
    X = model.matrix(~X, tmp)
    y = tmp$Norm_Value
    
    ### create design matrix
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
    
    ### fit model and get the summary
    fit = lm(y ~ X + 0)
    res = summary(fit)
    
    ### reduce the memory size
    res = stripGlmLR(res)
    return(res)
})

###################################################
# Save the results
###################################################

fdiry = file.path(FD_RES, "model_linear")
fname = paste0("res_interactive_", CHROM, ".rds") # fname="res_interactive_chr21.rds"
fpath = file.path(fdiry, fname)
saveRDS(lst_res, fpath)
print("Done")


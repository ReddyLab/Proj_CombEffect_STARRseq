### set env
source("/home/mount/project/config_sing.R")
#print(dir(FD_WORK))
#print(dir(FD_RES))

### import
ctypes = c(col_character(), col_integer(), col_integer(), col_integer(),
           col_character(), col_integer(), col_integer(),
           col_character(), col_double(),  col_integer())
cnames = c("Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag",
           "Chrom_MTF",  "Start_MTF",  "End_MTF",
           "Motif", "Score", "Overlap")

fdiry = file.path(FD_RES, "annotation_fragment", "Input1")
fname = "chr17.bed"
fpath = file.path(fdiry, fname)

dat = read_tsv(fpath, col_types=ctypes, col_names=cnames)
print(head(dat, 3))


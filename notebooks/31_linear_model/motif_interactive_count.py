###################################################
# Set environment
###################################################

### set environment
import argparse
import sys
sys.path.append("/home/mount/project")

###################################################
# Import packages and parse arguments
###################################################

### import packages and paths
from config_sing import *
import multiprocessing
from functools import reduce

### parse argument from command line
parser = argparse.ArgumentParser(description='Get count table from annotated fragment.')

parser.add_argument('-t', '--target',  
                    type = str,  
                    help = 'Target region')

parser.add_argument("-c", "--core", 
                    type    = int,
                    default = 1,
                    help    = "Number of cores (Default: 1)")

parser.add_argument("-o", "--outdir", 
                    type = str,
                    help = "Output file directory")
                    
parser.add_argument("-s", "--score", 
                    type    = float,
                    default = 0,
                    help    = "Threshold of motif score")

parser.add_argument("-v", "--verbose", 
                    action = "store_true",
                    help   = "Set output verbosity")

###################################################
# Set global variables
###################################################

### global variables from argument
args = parser.parse_args()
TARGET          = args.target
FD_OUT          = os.path.join(FD_RES, 'model_linear', args.outdir, args.target)
N_CORE          = args.core
THRESHOLD_MOTIF = args.score
VERBOSE         = args.verbose

### create output directory
os.makedirs(FD_OUT, exist_ok=True)

### samples
SAMPLES = np.concatenate([
    reduce(np.char.add, ["Input", np.arange(1,6).astype(str)        ]),
    reduce(np.char.add, ["Input", np.arange(1,6).astype(str), "_20x"]),
    reduce(np.char.add, ["TFX",   np.arange(2,6).astype(str), "_DMSO"]),
    reduce(np.char.add, ["TFX",   np.arange(2,6).astype(str), "_Dex"])
])

### all motif file names
sample = "Input1_20x"
fdiry  = os.path.join(FD_RES, "annotation_fragment", sample, TARGET)
fname  = "*_merge.bed.gz"
fpath  = os.path.join(fdiry, fname)
MOTIFS = np.sort([os.path.basename(fp) for fp in glob(fpath)])

### PRINT
print("Target:                  ", TARGET)
print("Outdir:                  ", FD_OUT)
print("N_Core:                  ", N_CORE)
print("Threshold of motif score:", THRESHOLD_MOTIF)
print("Is Verbose:              ", VERBOSE)

###################################################
# Helper function to get motif count table
###################################################

def get_motif_data(motif_pair, sample, threshold_score, verbose=True):
    """..."""
    ### INIT
    lst_dat = []
    is_miss = False
    mtfs    = [motif.replace("_merge.bed.gz", "") for motif in motif_pair]
    txt     = "_".join(mtfs)
    
    for mtf, motif in zip(mtfs, motif_pair):
        ### INIT
        lst   = []
        idx   = -1 # HANDLE EXCEPTION: Empty file
        
        ### set input file path
        fdiry = os.path.join(FD_RES, "annotation_fragment", sample, TARGET)
        fname = motif
        fpath = os.path.join(fdiry, fname)
        
        ### PRINT: results of filtering
        if verbose:
            print(txt, "Import", sample, mtf, fpath)
        
        ### import data
        with gzip.open(fpath, 'rb') as file:
            for idx, line in enumerate(file):

                ### preprocess each line 
                line = str(line, 'utf-8')
                line = line.strip().split("\t")

                ### extract needed values
                idx1, idx2, idx3, idx4 = 5, 6, 8, 9
                len_mtf = int(line[idx2]) - int(line[idx1])
                len_lap = int(line[idx4])
                score   = float(line[idx3])

                ### filtering motifs scores and make sure
                ### fragments cover the full motif
                if (len_mtf <= len_lap) & (score >= threshold_score):
                    lst.append(line)

        ### PRINT: results of filtering
        if verbose:
            print(txt, "Filter", sample, mtf, str(idx + 1) + "-" + str(len(lst)), flush=True)        
    
        ### HANDLE EXCEPTION: Empty data
        if len(lst) == 0:
            if verbose:
                print(txt, "Skip Empty", sample, mtf, flush=True) 
            is_miss = True
            break
        
        ### update list
        lst_dat.append(lst)
        
    
    ### concatenate motif names
    return lst_dat, txt, is_miss

def get_motif_table(lst, sample, verbose=True):
    """..."""
    ### wrap up a list of lines into a dataframe
    CNAMES = ["Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag",
              "Chrom_MTF",  "Start_MTF",  "End_MTF",
              "Motif",      "Score",      "Overlap"]
    CTYPES = [np.str, np.int,   np.int, np.int,
              np.str, np.int,   np.int,
              np.str, np.float, np.int]
    dat = pd.DataFrame(lst, columns=CNAMES) \
        .astype(dict(zip(CNAMES, CTYPES))) \
        .assign(Sample = sample)

    ### select columns
    dat = dat.loc[:,["Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag", "Sample", "Motif"]]
    return dat

def fun_motif(motif1, motif2):
    """..."""
    motif = motif1.astype(str) + "_" + motif2.astype(str)
    motif = motif.str.replace("_nan|nan_", "")
    return motif

def fun_design(sample, motif):
    """..."""
    group = sample.str.replace("[0-9]|_20x", "")
    if "Input" == group[0]:
        return "Input"
    else:
        return group.astype(str) + "_" + motif.astype(str)

def get_motif_count(motif_pair, verbose=True):
    """..."""
    ### INIT
    is_created = False
    
    for sample in SAMPLES:
        ###################################################
        ### Get motif table
        ###################################################
        
        ### import data for each motif, skip if any is empty
        lst_dat, txt, is_miss = get_motif_data(motif_pair, sample, THRESHOLD_MOTIF, verbose=verbose)
        if is_miss:
            continue

        ### arrange data into dataframe
        df1 = get_motif_table(lst_dat[0], sample, verbose=verbose)
        df2 = get_motif_table(lst_dat[1], sample, verbose=verbose)
        lst_dat = None
        
        ###################################################
        ### Get motif count
        ###################################################
        
        ### merge motif annotation
        dat = pd.merge(df1, df2, on=['Chrom_Frag', 'Start_Frag', "End_Frag", "Count_Frag", "Sample"], how='outer')
        dat = dat.assign(Motif = lambda x: fun_motif(x.Motif_x, x.Motif_y))
        dat = dat.assign(X     = lambda x: fun_design(x.Sample, x.Motif))

        ### summarize into counts for each motif combination
        dat = dat.astype({"Motif_x": str, "Motif_y": str, "Motif": str})
        dat = dat  \
            .groupby(["Sample", "Motif_x", "Motif_y", "X"]) \
            .agg(Value=('Count_Frag', sum)) \
            .reset_index()
            
        ###################################################
        ### Store the count table
        ###################################################
        
        ### SET: file path for output count table
        fdiry = FD_OUT 
        fname = "count_" + txt + ".tsv"
        fpath = os.path.join(fdiry, fname)
        
        ### store results
        ### create the table for the first sample or if the file is not yet created
        ### Otherwise, append the counts in the file
        if is_created:

            ### PRINT: file path for output count
            if verbose:
                print(txt, sample, "Store_Append", fpath, flush=True)

            ### append the file
            dat.to_csv(fpath, sep='\t', index=False, mode='a', header=False)
            
        else:

            ### PRINT: file path for output count
            if verbose:
                print(txt, sample, "Store_Create", fpath, flush=True)
            
            ### create the file
            dat.to_csv(fpath, sep='\t', index=False, mode='w')

            ### update flag
            is_created = True

###################################################
# Get motif count table
###################################################

### Test
#motif_pair = ("AP1_1_merge.bed.gz", "AP1_2_merge.bed.gz")
#get_motif_count(motif_pair, verbose=VERBOSE)
#motif_pairs = it.combinations(MOTIFS[:3], 2)

### INIT
motif_pairs = it.combinations(MOTIFS, 2)

### run analysis
pool = multiprocessing.Pool(N_CORE)
res  = pool.map(get_motif_count, motif_pairs)



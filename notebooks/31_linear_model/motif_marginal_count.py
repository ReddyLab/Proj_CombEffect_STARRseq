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
sam    = "Input1_20x"
fdiry  = os.path.join(FD_RES, "annotation_fragment", sam, TARGET)
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
# Set global variables
###################################################

def get_count_table(motif, verbose = True):
    """import annotated fragment and generate count table"""
    
    ### INIT
    mtf = motif.replace("_merge.bed.gz", "")
    is_created = False

    ### for each sample, import annotated fragment and generate count table
    for sam in SAMPLES:    
        
        ###################################################
        ### Import & filtering fragment annotation
        ###################################################
        
        ### set input file path
        fdiry = os.path.join(FD_RES, "annotation_fragment", sam, TARGET)
        fname = motif
        fpath = os.path.join(fdiry, fname)
        if verbose:
            print(mtf, sam, "Import", fpath, flush=True)
        
        ### INIT
        lst = []
        idx = -1 # HANDLE EXCEPTION: Empty file
        
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
                if (len_mtf <= len_lap) & (score >= THRESHOLD_MOTIF):
                    lst.append(line)
        
        if verbose:
            print(mtf, sam, "Filter", str(idx + 1) + "-" + str(len(lst)), flush=True) 
            
        ### HANDLE EXCEPTION: Empty data
        if len(lst) == 0:
            if verbose:
                print(mtf, sam, "Skip Empty", flush=True) 
            continue
        
        ###################################################
        ### summarize the annotation and create count table
        ###################################################
        
        ### wrap up a list of lines into a dataframe
        CNAMES = ["Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag",
                  "Chrom_MTF",  "Start_MTF",  "End_MTF",
                  "Motif",      "Score",      "Overlap"]
        CTYPES = [np.str, np.int,   np.int, np.int,
                  np.str, np.int,   np.int,
                  np.str, np.float, np.int]
        dat = pd.DataFrame(lst, columns=CNAMES) \
            .astype(dict(zip(CNAMES, CTYPES))) \
            .assign(Sample = sam)
        
        ### count the number of motifs for each fragment
        dat = dat  \
            .groupby(["Chrom_Frag", "Start_Frag", "End_Frag", "Count_Frag", "Motif", "Sample"]) \
            .size() \
            .reset_index(name='N_Motif')
        
        ### Summarize into count table
        dat = dat \
            .groupby(["Sample", "Motif", "N_Motif"]) \
            .agg(Value=('Count_Frag', sum)) \
            .reset_index()
        
        if verbose:
            print(mtf, sam, "Count", np.sum(dat.Value), flush=True) 
        
        ###################################################
        ### Store the count table
        ###################################################
        
        ### SET: file path for output count table
        fdiry = FD_OUT 
        fname = "count_" + mtf + ".tsv"
        fpath = os.path.join(fdiry, fname)

        ### store results
        ### create the table for the first sample or if the file is not yet created
        ### Otherwise, append the counts in the file
        if is_created:

            ### PRINT: file path for output count
            print(mtf, sam, "Store_Append", fpath, flush=True)

            ### append the file
            dat.to_csv(fpath, sep='\t', index=False, mode='a', header=False)
            
        else:

            ### PRINT: file path for output count
            print(mtf, sam, "Store_Create", fpath, flush=True)
            
            ### create the file
            dat.to_csv(fpath, sep='\t', index=False, mode='w')

            ### update flag
            is_created = True
    

###################################################
# Get motif count table for marginal analysis
###################################################
print("\n++++++++++ Start Analysis ++++++++++")

### INIT: timer
tic = time.time()

### get countable table and store
pool = multiprocessing.Pool(N_CORE)
res  = pool.map(get_count_table, MOTIFS)

### PRINT: end message
toc = time.time()
print("Script Done!")
print("Time Elapse:" , str(timedelta(seconds = toc - tic)))


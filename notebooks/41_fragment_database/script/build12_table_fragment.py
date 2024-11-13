"""
The script is to create the table of fragment
"""

#######################################################
### Set environment
###++++++++++++++++++++++++++++++++++++++++++++++++++++

### import common packages
import numpy  as np
import itertools as it
import sys, os, gzip
from   functools import reduce

### update print
from functools import partial
print = partial(print, flush=True)

### set working directories
sys.path.append('/mount/project')
from config_sing_path import *

### import specific packages
import sqlite3
# https://stackoverflow.com/questions/49456158/integer-in-python-pandas-becomes-blob-binary-in-sqlite
sqlite3.register_adapter(np.int64, lambda val: int(val))
sqlite3.register_adapter(np.int32, lambda val: int(val))

#######################################################
### Parse arguments
###++++++++++++++++++++++++++++++++++++++++++++++++++++
import argparse

### set arguments
parser = argparse.ArgumentParser()
parser.add_argument("--chrom",   type = str, help = "Chromosome")
parser.add_argument("--prefix",  type = str, help = "Prefix of the database name", default="fragment")
parser.add_argument("--fout",    type = str, help = "Output directory")
parser.add_argument("--finp",    type = str, help = "Input directory")
parser.add_argument("--verbose", type = str, help = "Should the processing message be printed?", default=True)

### parse arguments
args    = parser.parse_args()
CHROM   = args.chrom
FD_OUT  = args.fout
FD_INP  = args.finp
PREFIX  = args.prefix
VERBOSE = args.verbose

#######################################################
### Global varialbes and I/O
###++++++++++++++++++++++++++++++++++++++++++++++++++++

### file path of fragment database
fdiry  = FD_OUT
fname  = f"{PREFIX}_{CHROM}.db"
FP_DTB = os.path.join(fdiry, fname)

### Set Samples
fun = np.core.defchararray.add
idx = np.arange(1,6).astype("str")

INPUT    = reduce(fun, ["Input", idx])
INPUT20X = reduce(fun, ["Input", idx,     "_20x"])
TFX_DMSO = reduce(fun, ["TFX",   idx[1:], "_DMSO"])
TFX_DEX  = reduce(fun, ["TFX",   idx[1:], "_Dex"])
SAMPLES  = np.concatenate([INPUT, INPUT20X, TFX_DMSO, TFX_DEX])

### show info
if (VERBOSE):
    print("Global variables:")
    print("Chromsome: ", CHROM)
    print("Database:  ", FP_DTB)
    print("Input file directory:", FD_INP)
    print()

##################################################
### Helper functions
### ++++++++++++++++++++++++++++++++++++++++++++++

### helper function to process each row
def prep_line(line):
    """Function to process each line"""
    ### Decode
    lst = line.decode('ASCII').strip().split('\t') 

    ### parse info
    key = "_".join(lst[0:3])
    val = lst[0:3] + lst[4:-1]
    return [key] + val

### helper function to get a chunk of file
def get_chunks(gen, rows=10000):
    """Divides the data into 10000 rows each """
    iterable = iter(gen)
    while True:
        x = list(it.islice(iterable, rows))
        if not x:
            return
        yield x

##################################################
### 
### ++++++++++++++++++++++++++++++++++++++++++++++

fdiry = os.path.join(FD_INP, "Input1_20x")
fname = "target_PER1.bed.gz"
fpath = os.path.join(fdiry, fname)

n_lines  = 5
n_chunks = 3

with gzip.open(fpath, "rb") as file:
    header = file.readline()
    chunks = get_chunks(file, rows=n_lines)
    chunks = it.islice(chunks, n_chunks)

    for chunk in chunks:
        
##################################################
### Set environment
### ++++++++++++++++++++++++++++++++++++++++++++++

### basic
import sys
sys.path.append('/mount/project')
from config_sing import *

### update print
from functools import partial
print = partial(print, flush=True)

### specific tools
from functools import reduce
import itertools as it
import sqlite3
# https://stackoverflow.com/questions/49456158/integer-in-python-pandas-becomes-blob-binary-in-sqlite
sqlite3.register_adapter(np.int64, lambda val: int(val))
sqlite3.register_adapter(np.int32, lambda val: int(val))

### parse argument
import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    'chrom', 
    type=str, 
    help='Chromsome')
args = parser.parse_args()

##################################################
### Global variables
### ++++++++++++++++++++++++++++++++++++++++++++++

### chromsome
CHROM = args.chrom

### file path of fragment database
fdiry = os.path.join(FD_RES, 'database')
fname = f"fragment_{CHROM}.db"
FPATH_DB = os.path.join(fdiry, fname)

### set name of table files
FNAME_TXT = f"{CHROM}.bed.gz"

### set Samples
fun = np.core.defchararray.add
idx = np.arange(1,6).astype("str")
INPUT    = reduce(fun, ["Input", idx])
INPUT20X = reduce(fun, ["Input", idx,     "_20x"])
TFX_DMSO = reduce(fun, ["TFX",   idx[1:], "_DMSO"])
TFX_DEX  = reduce(fun, ["TFX",   idx[1:], "_Dex"])
SAMPLES  = np.concatenate([INPUT, INPUT20X, TFX_DMSO, TFX_DEX])

### show info
print("Global variables:")
print("Chromsome:      ", CHROM)
print("Database:       ", FPATH_DB)
print("Table file name:", FNAME_TXT)
print()

##################################################
### Helper function
### ++++++++++++++++++++++++++++++++++++++++++++++

### helper function to process each row
def prep_line(lst):
    """..."""
    ### parse info and get fragment and motif ID
    fragment = "_".join(lst[:3])
    binding  = "_".join(lst[4:8])
    overlap  = int(lst[-1])
    
    ### calc motif lengths
    mtf_len  = int(lst[6]) - int(lst[5])
    
    return mtf_len == overlap, (fragment, binding)

##################################################
### Set SQL query
### ++++++++++++++++++++++++++++++++++++++++++++++

query_reset_table = "DROP TABLE IF EXISTS Annotation"
query_reset_index = "DROP INDEX IF EXISTS idx_annot_frag"

query_table = """
    CREATE TABLE IF NOT EXISTS Annotation (
        fragment TEXT, 
        binding  TEXT,
        FOREIGN KEY (fragment) REFERENCES Fragment (fragment),
        FOREIGN KEY (binding)  REFERENCES Motif    (binding),
        UNIQUE (fragment, binding) ON CONFLICT IGNORE
    );"""

query_index = "CREATE INDEX idx_annot_frag ON Annotation (fragment)"

query_insert = """
    INSERT OR IGNORE INTO Annotation
        (fragment, binding)
    VALUES 
        (?,?)
    """

##################################################
### Create Annotation Table
### ++++++++++++++++++++++++++++++++++++++++++++++

### set data directory
fdirys = list()
samples = "Input?,Input?_20x,TFX?_DMSO,TFX?_Dex".split(",")
for sam in samples:
    fglob = os.path.join(FD_RES, "annotation", sam)
    lst   = glob.glob(fglob)
    lst.sort()
    fdirys += lst
fname = FNAME_TXT
fpaths_gz = [os.path.join(fdiry, fname) for fdiry in fdirys]
fpath_db  = FPATH_DB

### create table
counter_tot = 0
with sqlite3.connect(fpath_db) as conn:
    ### reset table
    cursor = conn.cursor()
    query  = query_reset_table
    cursor.execute(query) 
    query  = query_reset_index
    cursor.execute(query)
    
    ### create table
    cursor = conn.cursor()
    query  = query_table
    cursor.execute(query)
    query  = query_index
    cursor.execute(query)
    
    ### loop through each annotation file
    for fpath_gz in fpaths_gz:
        ### init, show progress
        print(fpath_gz)
        counter = 0
        
        with gzip.open(fpath_gz, "rb") as file:
            
            ### insert values
            query  = query_insert
            #n_line = 5
            #lines  = it.islice(file, n_line)
            lines  = file

            for idx, line in enumerate(lines):
                lst = line.decode('ASCII').strip().split('\t')  
                is_cover, row = prep_line(lst)
                counter_tot += is_cover
                counter     += is_cover
                if is_cover:
                    cursor.execute(query, row)
        print(f"#Rows: {idx+1}; Inserted: {counter}")
              
    print("\n#Rows Total:", counter_tot) 

    ### show created table info
    cursor = conn.cursor()
    query = "select count(*) from Annotation"
    cursor.execute(query)
    print("#Rows Table:", cursor.fetchall())
    
    ### show that the table is created
    #cursor = conn.cursor()
    #cursor.execute("SELECT * FROM Annotation")
    #for row in cursor.fetchall():
    #    print(row)
    
    
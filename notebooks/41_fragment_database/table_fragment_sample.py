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
def prep_line(line, sample):
    """..."""
    ### Decode
    lst = line.decode('ASCII').strip().split('\t') 

    ### parse info and get fragment and motif ID
    fragment = "_".join(lst[:3])
    count    = lst[-1]
    
    return fragment, sample, count

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
### Set SQL query
### ++++++++++++++++++++++++++++++++++++++++++++++

query_reset_table = "DROP TABLE IF EXISTS Count"
query_reset_index = "DROP INDEX IF EXISTS idx_count_sample"

query_table = """
    CREATE TABLE IF NOT EXISTS Count (
        fragment TEXT, 
        sample   TEXT,
        count    INTEGER,
        FOREIGN KEY (fragment) REFERENCES Fragment (fragment),
        FOREIGN KEY (sample)   REFERENCES Sample   (sample)
    );"""

query_index = "CREATE INDEX idx_count_sample ON Count (sample)"

query_insert = ("""INSERT OR IGNORE INTO Count
    (fragment, sample, count)
    VALUES 
    (?,?,?)""")

##################################################
### Create Count Table
### ++++++++++++++++++++++++++++++++++++++++++++++

### set data directory
fdirys_gz = list()
samples = "Input?,Input?_20x,TFX?_DMSO,TFX?_Dex".split(",")
for sam in samples:
    fglob = os.path.join(FD_RES, "count_fragment", sam)
    lst   = glob.glob(fglob)
    lst.sort()
    fdirys_gz += lst
fname_gz = FNAME_TXT

fpath_db = FPATH_DB

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
    for fdiry_gz in fdirys_gz:
        
        ### init
        sam = os.path.basename(fdiry_gz)
        fpath_gz = os.path.join(fdiry_gz, fname_gz)
        counter = 0

        ### show progress
        print(fpath_gz)
        print("Sample: ", sam)
        
        
        with gzip.open(fpath_gz, "rb") as file:
            
            ### set query
            query  = query_insert

            ### set lines and chunks
            #n_lines = 10
            #lines  = it.islice(file, n_lines)
            lines   = file
            n_chunk = 1000000
            chunks  = get_chunks(lines, rows=n_chunk)

            for chunk in chunks:
                ### preprocess lines
                fun = partial(prep_line, sample=sam)
                lst = list(map(fun, chunk))

                ### insert
                cursor.executemany(query, lst)

                ### count
                counter_tot += len(lst)
                counter     += len(lst)

            ### show progress for each sample                
            print(f"#Rows: {counter}")

    ### show total number of lines read from files
    print("\n#Rows Total:", counter_tot)

    ### show created table info
    cursor = conn.cursor()
    query = "select count(*) from Count"
    cursor.execute(query)
    print("#Rows Table:", cursor.fetchall())
    
    ### show that the table is created
    #cursor = conn.cursor()
    #cursor.execute("SELECT * FROM Count")
    #for row in cursor.fetchall():
    #    print(row)
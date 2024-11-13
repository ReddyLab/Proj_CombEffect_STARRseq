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
FD_RES = "/gpfs/fs1/data/reddylab/Kuei/out/proj_combeffect"

### import specific packages
import sqlite3
# https://stackoverflow.com/questions/49456158/integer-in-python-pandas-becomes-blob-binary-in-sqlite
sqlite3.register_adapter(np.int64, lambda val: int(val))
sqlite3.register_adapter(np.int32, lambda val: int(val))


### import packages for benchmark performance
import cProfile, pstats, time, timeit
import matplotlib.pyplot as plt

##################################################
### Helper functions
### ++++++++++++++++++++++++++++++++++++++++++++++

### helper function to get a chunk of file
def get_chunks(gen, rows=10):
    """Divides the data into chunks with size as rows"""
    iterable = iter(gen)
    while True:
        x = list(it.islice(iterable, rows))
        if not x:
            return
        yield x

### helper function to process each row
def prep_line(line):
    """Function to process each line"""
    ### Decode
    lst = line.decode('ASCII').strip().split('\t') 

    ### parse info
    key = "_".join(lst[0:3])
    val = lst[0:3] + lst[4:-1]
    return [key] + val

def gen_line(file, n_chunksize=None, n_lines=None):
    """generate lines or chunks of lines from the file"""
    ### remove file header
    header = file.readline()
    lines  = file
    
    ### preprocess each line
    fun = prep_line
    gen = map(fun, lines)
    
    ### set number of lines read if specified
    if n_lines is not None:
        gen = it.islice(gen, n_lines)
    
    ### set chunks if specified
    if n_chunksize is not None:
        gen = get_chunks(gen, n_chunksize)

    return gen

##################################################
### Set SQL query
### ++++++++++++++++++++++++++++++++++++++++++++++

query_reset_table = "DROP TABLE IF EXISTS Fragment"

query_table_frag = ("""
    CREATE TABLE IF NOT EXISTS Fragment(
        fragment TEXT PRIMARY KEY, 
        chrom    TEXT,
        start    INTEGER,
        end      INTEGER,
        pct_at   REAL,
        pct_gc   REAL,
        num_A    INTEGER,
        num_C    INTEGER,
        num_G    INTEGER,
        num_T    INTEGER,
        num_N    INTEGER,
        num_oth  INTEGER
    );""")

query_table_auto = ("""
    CREATE TABLE IF NOT EXISTS Fragment(
        fragment TEXT, 
        chrom    TEXT,
        start    INTEGER,
        end      INTEGER,
        pct_at   REAL,
        pct_gc   REAL,
        num_A    INTEGER,
        num_C    INTEGER,
        num_G    INTEGER,
        num_T    INTEGER,
        num_N    INTEGER,
        num_oth  INTEGER
    );""")

query_insert = ("""
    INSERT OR IGNORE INTO Fragment
        (fragment, chrom, start, end, pct_at, pct_gc,
         num_A, num_C, num_G, num_T, num_N, num_oth) 
    VALUES 
        (?,?,?,?,?,?,?,?,?,?,?,?)
    """)

##################################################
### Set database function
### ++++++++++++++++++++++++++++++++++++++++++++++

def refresh(query_table, fpath_database):
    """
    Helper function to refresh the database by 
    deleting original table and create a new one
    """
    with sqlite3.connect(fpath_database) as conn:
        ### init
        cursor = conn.cursor()

        ### reset table
        query  = query_reset_table
        cursor = cursor.execute(query)

        ### create table
        query  = query_table
        cursor = cursor.execute(query)
                
##################################################
### Helper function to check if the table is created correctly
### ++++++++++++++++++++++++++++++++++++++++++++++         
                
def check_table_size(n_lines, fpath_database):
    """count the number of rows/lines in a table created in a database"""
    with sqlite3.connect(fpath_database) as conn:
        ### initiation
        cursor = conn.cursor()
        query  = "select count(*) from Fragment"
        
        ### get the table size of table
        cursor = cursor.execute(query)
        counts = cursor.fetchall()
        counts = counts[0][0]
    
    if counts == n_lines:
        print("Check table size:  passed!")    
    else:
        print("Check table size:  failed.")
    

def check_table_lines(fpath_database, fpath_table):
    """check the integrity of the table creation"""
    
    def get_line_from_database():
        """generator of lines from database"""
        with sqlite3.connect(fpath_database) as conn:
            ### initiation
            cursor = conn.cursor()
            query  = "select * from Fragment"

            ### get the table size of table
            cursor = cursor.execute(query)
            for line in cursor:
                yield line
                
    def get_line_from_file():    
        """generator of lines from file"""
        with gzip.open(fpath_table, "rb") as file:
            lines = gen_line(file)
            for line in lines:
                yield line
    
    ### compare line by line
    lines_base = get_line_from_database()
    lines_file = get_line_from_file()
    fun = lambda x, y: str(x) == str(y)
    
    for line_base, line_file in zip(lines_base, lines_file):
        ### compare elements for each pair of lines and 
        ### return a list of true/false whether each pair of elements are equal
        res = list(it.starmap(fun, zip(line_base, line_file)))
        
        ### if all elements are equal for each line, continue
        ### if not, exit the function and print the failed message
        if all(res):
            continue
        else:
            print("Check table lines: failed.")
            return
        
    ### all lines are equal
    print("Check table lines: passed!")    
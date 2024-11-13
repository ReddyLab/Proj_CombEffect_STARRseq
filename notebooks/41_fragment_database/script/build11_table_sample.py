"""
The script is to create the table of sample
"""

#######################################################
### Set environment
###++++++++++++++++++++++++++++++++++++++++++++++++++++

### import common packages
import numpy  as np
import pandas as pd
import itertools as it
import sys, os
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
parser.add_argument("--fout",    type = str, help = "Output directory")
parser.add_argument("--finp",    type = str, help = "Input file of sample information and library size")
parser.add_argument("--prefix",  type = str, help = "Prefix of the database name", default="fragment")
parser.add_argument("--verbose", type = str, help = "Should the processing message be printed?", default=True)

### parse arguments
args    = parser.parse_args()
CHROM   = args.chrom
FD_OUT  = args.fout
FP_INP  = args.finp
PREFIX  = args.prefix
VERBOSE = args.verbose

#######################################################
### Global varialbes and I/O
###++++++++++++++++++++++++++++++++++++++++++++++++++++

### file path of fragment database
fdiry  = FD_OUT
fname  = f"{PREFIX}_{CHROM}.db"
FP_DTB = os.path.join(fdiry, fname)

### show info
if (VERBOSE):
    print("Global variables:")
    print("Chromsome:  ", CHROM)
    print("Database:   ", FP_DTB)
    print("Sample file:", FP_INP)
    print()

#######################################################
### Import data
###++++++++++++++++++++++++++++++++++++++++++++++++++++

dat_sam = pd.read_table(FP_INP)

if (VERBOSE):
    print("Import sample file:")
    print(dat_sam)
    print()

#######################################################
### Create table schema and queries
###++++++++++++++++++++++++++++++++++++++++++++++++++++

query_reset = ("DROP TABLE IF EXISTS Sample")

query_table = ("""CREATE TABLE IF NOT EXISTS Sample(
    sample    TEXT PRIMARY KEY, 
    treatment TEXT,
    size      INTEGER
);""")

query_insert = ("INSERT INTO Sample (sample, treatment, size) VALUES (?, ?, ?)")

#######################################################
### Create table 
###++++++++++++++++++++++++++++++++++++++++++++++++++++

with sqlite3.connect(FP_DTB) as conn:
    ### reset table
    cursor = conn.cursor()
    query  = query_reset
    cursor = cursor.execute(query)
    
    ### create table
    query  = query_table
    cursor = cursor.execute(query)
    
    ### insert values
    rows   = dat_sam.to_records(index=False)
    query  = query_insert
    cursor.executemany(query, rows)

#######################################################
### Check table
###++++++++++++++++++++++++++++++++++++++++++++++++++++

if (VERBOSE):
    print("Check if the table is created:")
    print("Show the table:")
    with sqlite3.connect(FP_DTB) as conn:
        ### show that the table is created
        cursor = conn.cursor()
        cursor = cursor.execute("SELECT * FROM Sample")
        for row in cursor.fetchall():
            print(row)
    print()
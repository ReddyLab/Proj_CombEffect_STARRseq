### basic
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, time, gzip, glob

#import itertools as it
#import os, sys, time
#from glob     import glob
#from datetime import timedelta
#from collections import Counter

###
SERVER = "Singularity | singularity_proj_combeffect"
#SERVER = "Duke Server: DCC"
#SERVER = "Duke Server: HARDAC"

### set paths
FD_WORK = "/mount/work"
FD_RLAB = "/mount/reddylab"
FD_PRJ  = "/mount/project"

FD_SRC  = os.path.join(FD_WORK, "source")
FD_EXE  = os.path.join(FD_WORK, "exe")
FD_ANN  = os.path.join(FD_WORK, "annotation")
FD_RES  = os.path.join(FD_WORK, "out", "proj_combeffect")

def show_env():
    print(f"You are in: {SERVER}")
    print(f"    BASE DIRECTORY:     {FD_WORK}") 
    print(f"    PATH OF SOURCE:     {FD_SRC}")
    print(f"    PATH OF EXECUTABLE: {FD_EXE}")
    print(f"    PATH OF ANNOTATION: {FD_ANN}")
    print(f"    PATH OF PROJECT:    {FD_PRJ}")
    print(f"    PATH OF RESULTS:    {FD_RES}")
    print()
    print("Library imported:")
    print("    numpy, pandas, matplotlib.pyplot")
    print("    os, sys, time, gzip, glob")
    print()
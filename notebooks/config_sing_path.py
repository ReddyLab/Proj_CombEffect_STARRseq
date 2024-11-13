"""
Setup the directories when running codes in singularity images of this project
This is a script to ensure that all scripts are using the same directories
"""

import os

### set main directories
FD_WORK = "/mount/work"
FD_RLAB = "/mount/reddylab"
FD_PRJ  = "/mount/project"

### set relative working directories
FD_SRC  = os.path.join(FD_WORK, "source")
FD_EXE  = os.path.join(FD_WORK, "exe")
FD_ANN  = os.path.join(FD_WORK, "annotation")
FD_RES  = os.path.join(FD_WORK, "out", "proj_combeffect")

### helpder function to check directories
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
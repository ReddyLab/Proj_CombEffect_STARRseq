### basic
import numpy     as np
import pandas    as pd
import os, sys, gzip

#import itertools as it
#import os, sys, time
#from glob     import glob
#from datetime import timedelta
#from collections import Counter

### set apths
FD_WORK = "/home/mount/work"
FD_RLAB = "/home/mount/reddylab"
FD_PRJ  = "/home/mount/project"

FD_SRC  = os.path.join(FD_WORK, "source")
FD_EXE  = os.path.join(FD_WORK, "exe")
FD_ANN  = os.path.join(FD_WORK, "annotation")
FD_RES  = os.path.join(FD_WORK, "out", "proj_combeffect")
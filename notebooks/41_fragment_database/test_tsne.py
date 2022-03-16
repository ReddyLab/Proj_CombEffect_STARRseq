### load tools
import os
from sklearn.datasets import make_blobs
from sklearn.manifold import TSNE
from matplotlib import pyplot as plt

###
import multiprocessing
print(multiprocessing.cpu_count())


### update print
from functools import partial
print = partial(print, flush=True)

### parse argument
import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    'n_core', 
    type=int, 
    help='Number of cores')
args = parser.parse_args()

N_CORE = args.n_cor

### generate data
n_samples  = 100000
n_centers  = 5
n_features = 10
X, y = make_blobs(
    n_samples    = n_samples, 
    n_features   = n_features, 
    cluster_std  = 1.0,
    centers      = n_centers, 
    shuffle      = False, 
    random_state = 42)

###
print("Generate Data:")
print(X.shape)
print(y.shape)

###
embeddings = TSNE(n_jobs=N_CORE, verbose=1).fit_transform(X)


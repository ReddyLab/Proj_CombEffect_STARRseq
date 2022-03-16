### load tools
import os
from sklearn.datasets import make_blobs
from MulticoreTSNE import MulticoreTSNE as TSNE
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


"""
fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(7, 3))

vis_x = X[:, 0]
vis_y = X[:, 1]
ax = axes[0]
im = ax.scatter(vis_x, vis_y, c=y, cmap=plt.cm.get_cmap("jet", n_centers), marker='.')

vis_x = embeddings[:, 0]
vis_y = embeddings[:, 1]
ax = axes[1]
im = ax.scatter(vis_x, vis_y, c=y, cmap=plt.cm.get_cmap("jet", n_centers), marker='.')

fig.subplots_adjust(right=0.8)
cax  = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cax)

fname = "mctsne_blob_" + str(N_CORE) + ".png"
plt.savefig(fname)
"""
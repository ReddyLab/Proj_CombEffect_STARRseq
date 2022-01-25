### set environment
import sys
sys.path.append('/home/mount/project/')
from config_sing import *

from functools import partial
print = partial(print, flush=True)

### import basic 
import tensorflow as tf
import json
from sklearn.model_selection import train_test_split

### from https://github.com/const-ae/Neural_Network_DNA_Demo
sys.path.append('Neural_Network_DNA_Demo/')
from helper import IOHelper, SequenceHelper

### import deepstarr
import DeepSTARR

### set GPU for tensorflow
physical_devices = tf.config.experimental.list_physical_devices("GPU")
tf.config.experimental.set_memory_growth(physical_devices[0], True)
print(physical_devices)

### get arguments
#import argparse
#parser=argparse.ArgumentParser()
#parser.add_argument('--params', help='file path of model parameters')
#parser.add_argument('--fout',   help='output directory')
#parser.add_argument('--prefix', help='prefix of file name')
#args=parser.parse_args()



print("load data")
fdiry = "/home/mount/work/out/proj_combeffect/peak/cradle_deepstarr_data"
fname = "whole_genome_X_sub.npy"
fpath = os.path.join(fdiry, fname)
with open(fpath, 'rb') as file:
    X = np.load(file)
print(X.shape)

print("load model")
fdiry = "/home/mount/work/out/proj_combeffect/peak"
model_name = "Model_DeepSTARR"
model_ID = os.path.join(fdiry, "cradle_deepstarr_results", model_name)
keras_model, keras_model_weights, keras_model_json = DeepSTARR.load_model(model_ID)
print(keras_model.summary())

print("\nRunning DeepExplain ...\n")
TASK1="DMSO"
TASK2="Dex"
class_output = TASK1
scores_task1 = DeepSTARR.my_deepExplainer(keras_model, X, class_output=class_output)

class_output = TASK2
scores_task2 = DeepSTARR.my_deepExplainer(keras_model, X, class_output=class_output)

print(len(scores_task1))
print(scores_task1[0].shape)
print(scores_task1[1].shape)
print(len(scores_task2))
print(scores_task2[0].shape)
print(scores_task2[1].shape)

print("\nSaving ...\n")

import h5py
import os


#if (os.path.isfile(sequence_set+"_"+model_ID_out+"_"+class_output+"_contribution_scores.h5")):
#    os.remove(str(sequence_set+"_"+model_ID_out+"_"+class_output+"_contribution_scores.h5"))
#f = h5py.File(sequence_set+"_"+model_ID_out+"_"+class_output+"_contribution_scores.h5")

fdiry = "/home/mount/work/out/proj_combeffect/peak"
fpath = os.path.join(
    fdiry, 
    "cradle_deepstarr_results", 
    "WG_" + model_name + "_contribution_scores.h5")

print(f"Score file: {fpath}")

if (os.path.isfile(fpath)):
    os.remove(str(fpath))
f = h5py.File(fpath)

g = f.create_group("hyp_contrib_scores")
# save the hypothetical contribution scores
#g.create_dataset(class_output, data=scores[0])
g.create_dataset(TASK1, data=scores_task1[0])
g.create_dataset(TASK2, data=scores_task2[0])
print("Done hyp scores for " + class_output)

g = f.create_group("contrib_scores")
# save the actual contribution scores
#g.create_dataset(class_output, data=scores[1])
g.create_dataset(TASK1, data=scores_task1[1])
g.create_dataset(TASK2, data=scores_task2[1])
print("Done contr scores for " + class_output)

f.close()
print("Scores saved")
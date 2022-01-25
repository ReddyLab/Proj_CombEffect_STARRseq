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


#import keras
#from tensorflow import keras
#import tensorflow.keras.layers as kl
#from tensorflow.keras.layers import Conv1D, MaxPooling1D
#from tensorflow.keras.layers import Dropout, Reshape, Dense, Activation, Flatten
#from tensorflow.keras.layers import BatchNormalization, InputLayer, Input
#from tensorflow.keras import models
#from tensorflow.keras.models import Sequential, Model
#from tensorflow.keras.optimizers import Adam
#from tensorflow.keras.callbacks import EarlyStopping, History, ModelCheckpoint

### set GPU for tensorflow
physical_devices = tf.config.experimental.list_physical_devices("GPU")
tf.config.experimental.set_memory_growth(physical_devices[0], True)
print(physical_devices)


### get arguments
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--params', help='file path of model parameters')
parser.add_argument('--fout',   help='output directory')
parser.add_argument('--prefix', help='prefix of file name')
args=parser.parse_args()

MODEL_NAME = args.prefix
FD_OUT     = args.fout

print("Arguements")
for key, val in vars(args).items():
    print(f"{key}: {val}")

### import model parameteres
#fpath = "./model/test.json"
fpath  = args.params
with open(fpath, 'r') as file:
    txt    = file.read()
    params = json.loads(txt)

print("Parameters")
for key, val in params.items():
    print(f"{key}: {val}")

print("Load data")

fdiry = "/home/mount/work/out/proj_combeffect/peak/cradle_deepstarr_data"

fname = "whole_genome_X.npy"
fpath = os.path.join(fdiry, fname)
with open(fpath, 'rb') as file:
    X = np.load(file)

fname = "whole_genome_Y.npy"
fpath = os.path.join(fdiry, fname)
with open(fpath, 'rb') as file:
    Y = np.load(file)

X_train, X_valid, y_train, y_valid = train_test_split(
    X, Y, 
    test_size=0.1, 
    random_state=123)

print(X.shape)
print(Y.shape)
print(X_train.shape)
print(X_valid.shape)
print(y_train.shape)
print(y_valid.shape)
print(DeepSTARR.DeepSTARR(params)[0].summary())

Y_train = [y_train[:,0], y_train[:,1]]
Y_valid = [y_valid[:,0], y_valid[:,1]]
main_model, main_params = DeepSTARR.DeepSTARR(params)
main_model, my_history  = DeepSTARR.train(main_model, X_train, Y_train, X_valid, Y_valid, main_params)

###################################################
# Save the model
#++++++++++++++++++++++++++++++++++++++++++++++++++

# init
model_name = MODEL_NAME
fdiry      = FD_OUT

# save model structure
model_json = main_model.to_json()
fname = model_name + '.json'
fpath = os.path.join(fdiry, fname)
with open(fpath, "w") as json_file:
    json_file.write(model_json)
    
# save model weights
fname = model_name + '.h5'
fpath = os.path.join(fdiry, fname)
main_model.save_weights(fpath)

print("Model Saved!")
import numpy as np
import math
import time
import datetime

import os
import argparse

import matplotlib.pyplot as plt

import tensorflow as tf

from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import LearningRateScheduler

from tensorflow.keras.layers import Dense, Conv1D, BatchNormalization, Activation
from tensorflow.keras.layers import Input, Cropping1D, Flatten
from tensorflow.keras.models import Model
from tensorflow.python.ops import math_ops

from tensorboard.plugins.hparams import api as hp

from sklearn.metrics import r2_score
from scipy.stats import pearsonr, spearmanr

from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.model_selection import train_test_split

from tqdm import tqdm


def lr_schedule(epoch, lr):
    #lr = 0.001
    if epoch == 5:
        lr *= 0.5
    elif epoch == 7:
        lr *= 0.5 ** 2
    elif epoch == 9:
        lr *= 0.5 ** 3
    elif epoch >= 11:
        lr *= 0.5 ** 4
    print('Learning rate: ', lr)
    return lr


def lr_schedule_const(epoch):
    lr = 0.001
    return lr


def check_argv():

    parser = argparse.ArgumentParser()
    parser.add_argument('--kernel_size0', default=11, required=True, type=int, help='kernel size for convolution')
    parser.add_argument('--kernel_size2', default=11, required=True, type=int, help='kernel size for convolution')
    parser.add_argument('--no_stacks', default=2, required=True, type=int, help='no RB stacks')
    parser.add_argument('--dilation1', default=4, required=True, type=int, help='dilation rate 1')
    parser.add_argument('--dilation2', default=4, required=True, type=int, help='dilation rate 2')
    parser.add_argument('--batch', default=32, required=True, type=int, help='batch size')
    parser.add_argument('--const_lr', default=False, action='store_true')


    return parser.parse_args()


def binarize_labels(labels):
    #A549: [0.25, 0.5, 0.95, 0.98] [0.00000000e+00 3.53096163e-05 4.08277148e-02 1.01717487e-01] 
    # ~100K, 20K
    #293T: [0.25, 0.5, 0.95, 0.98] [0. 0.00223609 0.04327499 0.08719538]
    # ~100K, 40K
    #293T: [0.25, 0.5, 0.98, 0.99] [0. 0.00223609 0.08719538 0.12098964]
    # ~40K, 20K

    q = 4.08277148e-02
    labels_ = []
    for x in labels:
        if x<=q:
            labels_.append([0])
        elif x>q:
            labels_.append([1])
    return np.array(labels_)


def RB_block(inputs,
             num_filters=32,
             kernel_size=11,
             strides=1,
             activation='relu',
             dilation_rate=1):
    """1D Convolution-Batch Normalization-Activation stack builder
    """
    conv = Conv1D(num_filters,
                  kernel_size=kernel_size,
                  strides=strides,
                  padding='same',
                  dilation_rate=dilation_rate)
    x = inputs
    for _ in range(2):
        x = BatchNormalization()(x)
        x = Activation(activation)(x)
        x = conv(x)
    return x


def model_CNN_full(input_shape, no_stacks, kernel_size, dil):
    """Model builder
    """
    inputs = Input(shape=input_shape)

    # initiate 
    x = Conv1D(32, kernel_size=1, strides=1, padding='same', dilation_rate=1)(inputs)
    # another Conv on x before splitting
    y = Conv1D(32, kernel_size=1, strides=1, padding='same', dilation_rate=1)(x)

    d = [1, dil[0], dil[1]] #dilation
    #d = [1, 4, 5]
    #kernels = [9, 11, 15]
    kernels = [kernel_size[0], 11, kernel_size[1]]
    for i in range(no_stacks):
        for block in range(4):
            x = RB_block(x, num_filters=32, kernel_size=kernels[i], strides=1, activation='relu', dilation_rate=d[i])
        if i!=(no_stacks-1):
            y = tf.keras.layers.add([Conv1D(32, kernel_size=1, strides=1, padding='same', dilation_rate=1)(x), y])

    x = Conv1D(32, kernel_size=1, strides=1, padding='same', dilation_rate=1)(x)
    # adding up with what was shortcut from the prev layers
    x = tf.keras.layers.add([x, y])

    x = Conv1D(1, kernel_size=1, strides=1, padding='same', dilation_rate=1)(x)
    x = Flatten()(x)
    outputs = Dense(1, activation='sigmoid')(x)

    model = Model(inputs=inputs, outputs=outputs)

    return model


region = ['chr2','chr4','chr6','chr8','chr10', 'chr11', 'chr12', 'chr13',
          'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
          'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


args = check_argv()
cell = 'A549'

datadir = os.path.join(f'data_new/{cell}/inputs_numpy/train_newMeth')

if args.const_lr:
    checkpoint_path = f"models_new/A549_cls1_k4me3/checkpoint/"
else:
    checkpoint_path = f"models_new/A549_cls1_k4me3/checkpoint/"
# Create a callback that saves the model's weights

cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath=checkpoint_path,
                                                 monitor='val_loss', mode='min',
                                                 verbose=1, save_best_only=False,
                                                 save_freq='epoch')

#log_dir = os.path.join(model_dir, "/fit_log")\
if args.const_lr:
    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=f"models_new/A549_cls1_k4me3/log", histogram_freq=1)
    #tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=f"models/A549_meta_cls_opt_batch{args.batch}_k9_11_15_dl1_4_5/log", histogram_freq=1)
else:
    #tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=f"models/A549_meta_cls_opt_001lr_batch{args.batch}_k9_11_15_dl1_4_5/log", histogram_freq=1)
    tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=f"models_new/A549_cls1_k4me3/log", histogram_freq=1)

es_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=2)

NE = 1
lr_scheduler = LearningRateScheduler(lr_schedule)

input_shape=(5, 250)

model_cnn = model_CNN_full(input_shape=input_shape, no_stacks=args.no_stacks, kernel_size=[args.kernel_size0,args.kernel_size2], dil=[args.dilation1, args.dilation2])

adam = Adam(learning_rate=0.001)

model_cnn.compile(loss=tf.keras.losses.BinaryCrossentropy(),
                  optimizer=adam,
                  metrics=tf.keras.metrics.Accuracy())
ep = 0

class_weight = {0: 0.34759133131208897, 1: 35.69278}

for _ in range(1):
    
    for chrn in region:

        print(f'processing {chrn}...')
        with open(os.path.join(datadir, chrn+'_seqs.npy'), 'rb') as f: 
            seqs = np.load(f)
        with open(os.path.join(datadir, chrn+'_me3.npy'), 'rb') as f: 
            me3 = np.load(f)
        with open(os.path.join(datadir, chrn+'_labels.npy'), 'rb') as f: 
            labels = np.load(f)

        labels_bin = binarize_labels(labels)

        metadata = []

        for s,m3 in tqdm(zip(seqs,me3)):
            x = np.concatenate((s,m3.reshape(1, -1)), axis=0)
            metadata.append(x)
            
        #metadata = np.array([x.T for x in metadata])
        metadata = np.array(metadata)

        # INPUT SHAPE IS EVERYTHING EXCEPT THE NUMBER OF SAMPLES
        #input_shape = metadata.shape[1:]
        print('x_train shape:', metadata.shape)
        print(metadata.shape[0], 'train samples')
        print('y_train shape:', labels.shape)

        # TRAINING

        if args.const_lr:
            history = model_cnn.fit(metadata, labels_bin, initial_epoch=ep, epochs=ep+NE,
                                    callbacks=[cp_callback, tensorboard_callback, es_callback], 
                                    validation_split=0.1, batch_size=args.batch,
                                    shuffle=True, class_weight=class_weight)
        else:
            history = model_cnn.fit(metadata, labels_bin, initial_epoch=ep, epochs=ep+NE,
                            callbacks=[lr_scheduler, cp_callback, tensorboard_callback, es_callback], 
                            validation_split=0.1, batch_size=args.batch,
                            shuffle=True, class_weight=class_weight)
       
        ep += NE

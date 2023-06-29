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
    parser.add_argument('--kernel_size1', default=11, required=False, type=int, help='kernel size for convolution')
    parser.add_argument('--kernel_size2', default=11, required=False, type=int, help='kernel size for convolution')
    parser.add_argument('--no_stacks', default=2, required=False, type=int, help='no RB stacks')
    parser.add_argument('--dilation1', default=4, required=False, type=int, help='dilation rate 1')
    parser.add_argument('--dilation2', default=4, required=False, type=int, help='dilation rate 2')
    parser.add_argument('--batch', default=32, required=False, type=int, help='batch size')
    parser.add_argument('--const_lr', default=False, action='store_true')
    parser.add_argument('--model_name', required=True, type=str, help='model name prefix')
    parser.add_argument('--inputs', required=True, type=str, help='folder with inputs split by chromosome')

    return parser.parse_args()


def binarize_labels(labels, q=4.08277148e-02):
    #provide a cutoff for label binarization

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

datadir = args.inputs

if not os.path.exists('data/models'):
    os.makedirs('data/models')

checkpoint_path = f"data/models/{args.model_name}/checkpoint/"
cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath=checkpoint_path,
                                                 monitor='val_loss', mode='min',
                                                 verbose=1, save_best_only=False,
                                                 save_freq='epoch')

tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=f"models/{args.model_name}/log", histogram_freq=1)
es_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=2)

NE = 1
lr_scheduler = LearningRateScheduler(lr_schedule)

input_shape=(4, 1000)

model_cnn = model_CNN_full(input_shape=input_shape, no_stacks=args.no_stacks, kernel_size=[args.kernel_size1,args.kernel_size2], dil=[args.dilation1, args.dilation2])

adam = Adam(learning_rate=0.001)

model_cnn.compile(loss=tf.keras.losses.BinaryCrossentropy(),
                  optimizer=adam,
                  metrics=tf.keras.metrics.Accuracy())
ep = 0

class_weight = {0: 0.34759133131208897, 1: 35.69278}

for _ in range(1):
    
    for chrn in region:

        chr_data_avail = os.listdir(datadir)
        chr_data_avail = [x.split('.')[0].split('_')[0] for x in chr_data_avail]
        if chrn not in chr_data_avail:
            continue

        print(f'processing {chrn}...')
        with open(os.path.join(datadir, chrn+'_seqs.npy'), 'rb') as f: 
            seqs = np.load(f)
        with open(os.path.join(datadir, chrn+'_labels.npy'), 'rb') as f: 
            labels = np.load(f)

        labels_bin = binarize_labels(labels, q=4.08277148e-02)

        # TRAINING

        if args.const_lr:
            history = model_cnn.fit(seqs, labels_bin, initial_epoch=ep, epochs=ep+NE,
                                    callbacks=[cp_callback, tensorboard_callback, es_callback], 
                                    validation_split=0.1, batch_size=args.batch,
                                    shuffle=True, class_weight=class_weight)
        else:
            history = model_cnn.fit(seqs, labels_bin, initial_epoch=ep, epochs=ep+NE,
                            callbacks=[lr_scheduler, cp_callback, tensorboard_callback, es_callback], 
                            validation_split=0.1, batch_size=args.batch,
                            shuffle=True, class_weight=class_weight)
       
        ep += NE

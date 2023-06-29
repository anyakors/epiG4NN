from tqdm import tqdm

import os
import argparse
import sys

import numpy as np
import pandas as pd

import tensorflow as tf

from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import LearningRateScheduler

from tensorflow.keras.layers import Dense, Conv1D, BatchNormalization, Activation
from tensorflow.keras.layers import Input, Cropping1D, Flatten
from tensorflow.keras.models import Model
from tensorflow.python.ops import math_ops


def complementary(let):
    # A-T, C-G
    let_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return let_dict(let)


def reverse_complement(seq):
    return "".join([complementary(let) for let in seq[::-1]])


def get_seq_from_txt(fasta_file):
    seq_dict = {}
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq_name = line.strip()[1:]
                seq_dict[seq_name] = ""
            else:
                seq_dict[seq_name] += line.strip()
    return seq_dict


def hot_encode_seq(let):
    if let == "A":
        return np.array([1, 0, 0, 0])
    elif let == "T":
        return np.array([0, 1, 0, 0])
    elif let == "C":
        return np.array([0, 0, 1, 0])
    elif let == "G":
        return np.array([0, 0, 0, 1])
    elif let == "N":
        return np.array([0, 0, 0, 0])


def transform_input(inputs_):
    inputs = []
    # hot-encode
    for i in tqdm(range(len(inputs_))):
        # hot-encode seq
        x = np.array([hot_encode_seq(let) for let in inputs_[i]]).T
        inputs.append(x)

    return np.array(inputs)


def check_argv():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument(
        "--seqs",
        required=True,
        type=str,
        help="input .fa file with 1000nt G4 sequences",
    )
    parser.add_argument(
        "--epi",
        type=str,
        help="input .csv file with 1000nt epigenetic arrays",
    )
    parser.add_argument(
        "--model",
        type=str,
        help="model to use: seq, h3k4me3, atac, h3k27ac, h3k4me1, h3k36me3, h3k9me3",
    )
    parser.add_argument(
        "--output", required=True, type=str, help="output filename prefix; will be saved to output folder"
    )

    return parser.parse_args()


args = check_argv()

seqs = get_seq_from_txt(args.seqs)

if args.model!='seq':
    epi = pd.read_csv(args.epi,  header=None, index_col=0)
    epi_dict = {i:np.array(epi.loc[i]) for i in epi.index}


metadata = []
ids_to_predict = []

if args.model=='seq':
    for i, s in tqdm(zip(seqs.keys(), seqs.values())):
        s_ = np.array([hot_encode_seq(let) for let in s]).T
        metadata.append(s_)
        ids_to_predict.append(i)
else:
    for i, s in tqdm(zip(seqs.keys(), seqs.values())):
        if i in epi_dict.keys():
            s_ = np.array([hot_encode_seq(let) for let in s]).T
            ac = epi_dict[i]
            x = np.concatenate((s_, ac.reshape(1, -1)), axis=0)
            metadata.append(x)
            ids_to_predict.append(i)

metadata = np.array(metadata)
print("Inputs shape:", np.shape(metadata))

models_to_eval = os.listdir("data/models/")

for model in models_to_eval:
    if args.model.lower() in model:
        input_shape = metadata.shape[1:]
        model_cnn = tf.keras.models.load_model(
            os.path.join("data/models/", model, "checkpoint")
        )
        y_pred = model_cnn.predict(metadata)
        
df_pred = pd.DataFrame({'id': ids_to_predict, 'Probability': y_pred.reshape(len(metadata))})

if not os.path.exists('output'):
    os.makedirs('output')

df_pred.to_csv(os.path.join("output", args.output+'.csv'), index=None)
print('Predictions saved to:', os.path.join("output", args.output+'.csv'))
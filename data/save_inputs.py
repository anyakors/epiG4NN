import pandas as pd
import numpy as np
import argparse

import os
import json
from tqdm import tqdm

from Bio import SeqIO

from sklearn.preprocessing import StandardScaler, MinMaxScaler


def complementary(let):
    # A-T, C-G
    if let == 'A':
        return 'T'
    if let == 'T':
        return 'A'
    if let == 'C':
        return 'G'
    if let == 'G':
        return 'C'
    if let == 'N':
        return 'N'


def make_dict_meth(me):
    print('original length:', len(me))
    me.drop_duplicates(inplace=True)
    print('length after deleting duplicates:', len(me))
    me_dict = {}
    for key,start,end,meth in zip(me[3],me[1],me[2],me[5]):
        if key not in me_dict.keys():
            me_dict[key] = []
        if int(end)>int(start):
            if float(meth)>=0:
                me_dict[key].extend([float(meth) for _ in range(int(end)-int(start))])
            else:
                me_dict[key].extend([0.0 for _ in range(int(end)-int(start))])
    print('PQS before filtering:', len(me_dict.keys()))
    for key in list(me_dict.keys()):
        if len(me_dict[key])!=1000:
            _ = me_dict.pop(key)
    print('PQS after filtering:', len(me_dict.keys()))
    
    return me_dict


def hot_encode_seq(let):
    if let == 'A':
        return np.array([1, 0, 0, 0])
    elif let == 'T':
        return np.array([0, 1, 0, 0])
    elif let == 'C':
        return np.array([0, 0, 1, 0])
    elif let == 'G':
        return np.array([0, 0, 0, 1])
    elif let == 'N':
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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hg19", default="hg19.fa", required=True, type=str, help="path to hg19.fa"
    )

    return parser.parse_args()


args = check_argv()

min_max_scaler_labels = MinMaxScaler(feature_range=(0,1))

scores = pd.read_csv('data/g4_scores_final_norm.bed', delimiter='\t', header=None)

print(max(scores[2]))
scores[2][scores[2] < 0] = 0
scores[2] = min_max_scaler_labels.fit_transform(np.array(scores[2]).reshape(-1, 1))
scores.head()

scores_dict = {}

for key,value in zip(scores[0], scores[2]):
    scores_dict[key] = {}
    scores_dict[key]['score'] = value if float(value)>=0 else 0.0

print(len(scores_dict.keys()))

print("Reading PQS file")

PQS = pd.read_csv('data/PQS_padded.bed', delimiter='\t', header=None)

region = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
         'chr14', 'chr15', 'chr16','chr17', 'chr18', 'chr19',
         'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


seq_dict = {}

print("Reading hg19 file")

with open(args.hg19, mode='r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        identifier = record.id
        description = record.description
        if identifier in region:
            sequence = record.seq
            seq_dict[identifier] = str(sequence)

for chrn in region:

    dir_atac = 'data/epi_bychr/'

    epi_chr_data_avail = os.listdir(dir_atac)
    epi_chr_data_avail = [x.split('.')[0] for x in epi_chr_data_avail]
    if chrn not in epi_chr_data_avail:
        continue

    print("Fetching sequences for PQS...")

    PQS_dict = {}

    for name,c,s,e,strand in zip(PQS[3], PQS[0], PQS[1], PQS[2], PQS[5]):
        if c==chrn:
            PQS_dict[name] = {}
            PQS_dict[name]['chr'] = c
            PQS_dict[name]['strand'] = strand
            PQS_dict[name]['start'] = int(s)
            PQS_dict[name]['end'] = int(e)

    for key in tqdm(PQS_dict.keys()):
        if PQS_dict[key]['chr']==chrn:
            seq = seq_dict[PQS_dict[key]['chr']][PQS_dict[key]['start'] : PQS_dict[key]['end']].upper()
            if PQS_dict[key]['strand']=='-':
                seq = ''.join([complementary(let) for let in seq])[::-1]
            PQS_dict[key]['sequence'] = seq

    print(f'...Processing {chrn}...')

    atac = pd.read_csv(os.path.join(dir_atac,chrn+'.bed'), delimiter='\t', header=None)

    atac_dict = make_dict_meth(atac)
    
    for pqs in tqdm(list(scores_dict.keys())):
        if pqs in atac_dict.keys():
            scores_dict[pqs]['id'] = pqs
            scores_dict[pqs]['atac'] = atac_dict[pqs]

    input_seqs = []
    input_atac = []
    labels = []

    for pqs in tqdm(scores_dict.keys()):
        if 'id' in scores_dict[pqs].keys() and pqs in PQS_dict.keys() and 'sequence' in PQS_dict[pqs].keys():
            input_seqs.append(PQS_dict[pqs]['sequence'])
            input_atac.append(np.array(scores_dict[pqs]['atac']))
            labels.append(scores_dict[pqs]['score'])
     
    input_seqs = transform_input(input_seqs)
    input_atac = np.array(input_atac)
    labels = np.array(labels)

    print(np.shape(input_seqs), np.shape(labels))

    if not os.path.exists('data/inputs_numpy'):
        os.makedirs('data/inputs_numpy')

    if chrn in ['chr1', 'chr3', 'chr5', 'chr7', 'chr9']:
        writedir = os.path.join('data/inputs_numpy/test_1000nt')
    else:
        writedir = os.path.join('data/inputs_numpy/train_1000nt')
    
    if not os.path.exists(writedir):
        os.makedirs(writedir)

    with open(os.path.join(writedir, f'{chrn}_seqs.npy'), 'wb') as f: 
        np.save(f, input_seqs)
    with open(os.path.join(writedir, f'{chrn}_epi.npy'), 'wb') as f: 
        np.save(f, input_atac)
    with open(os.path.join(writedir, f'{chrn}_labels.npy'), 'wb') as f: 
        np.save(f, labels)

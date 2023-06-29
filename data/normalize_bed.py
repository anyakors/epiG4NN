import pandas as pd
import numpy as np


df = pd.read_csv('data/g4_scores_final.bed', header=None, index_col=None, sep='\t')
df['score'] = df[2]/df[2].max()

with open('data/g4_scores_final_norm.bed', 'w') as f:
    for pqs, strand, score in zip(df[0], df[1], df['score']):
        f.write(f'{pqs}\t{strand}\t{score}\n')

df = pd.read_csv('data/epi_final.bed', header=None, index_col=None, sep='\t')
df['score'] = df[5]/df[5].max()

with open('data/epi_final_norm.bed', 'w') as f:
    for chrn, st, en, pqs, strand, score in zip(df[0], df[1], df[2], df[3], df[4], df['score']):
        f.write(f'{chrn}\t{st}\t{en}\t{pqs}\t{strand}\t{score}\n')
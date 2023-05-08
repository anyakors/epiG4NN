<p align="center"><img width=50% src="https://github.com/anyakors/epiG4NN/blob/main/media/logo.png"></p>

# epiG4NN: G-quadruplex prediction in live cells using epigenetic data

This repository contains the scripts for the **epiG4NN** data preparation pipeline and the code for epiG4NN model training. 

## Overview

**epiG4NN** model uses the potential G4 sequence padded to 1000 nucleotides matched with the respective epigenetic data of choice.
<p align="center"><img width=60% src="https://github.com/anyakors/epiG4NN/blob/main/media/scheme.jpeg"></p>

## Installation

Clone the **epiG4NN** repository to a local folder. Use the file `requirements.txt` to create a `conda` environment:

```
conda create --name epig4nn --file requirements.txt
conda activate epig4nn
```

## Data preparation pipeline

Download the human genome fasta file and supply it to the data preparation script. PQS sequences found in hg19 are supplementaed in the `data/` folder and are used by default.

```
cd data/
bash data_prep.sh
```

## Training

After generating the inputs with the `data_prep.sh` script, one can train either a bare sequence model (`train_seq.py`) or a model with epigenetic inputs (`train_epi.py`).

## Utils
`PQS_search.py` file contains a python re implementation for PQS search in the human genome

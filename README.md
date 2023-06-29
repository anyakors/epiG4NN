<p align="center"><img width=50% src="https://github.com/anyakors/epiG4NN/blob/main/media/logo.png"></p>

# epiG4NN: G-quadruplex prediction in live cells using epigenetic data

This repository contains the scripts for the **epiG4NN** data preparation pipeline and the code for epiG4NN model training. 

## Overview

**epiG4NN** model uses the potential G4 sequence padded to 1000 nucleotides matched with the respective epigenetic data of choice.
<p align="center"><img width=60% src="https://github.com/anyakors/epiG4NN/blob/main/media/scheme.jpeg"></p>

## Installation

Clone the **epiG4NN** repository to a local folder. Use the file `environment.yml` to create a `conda` environment:

```
conda env create -f environment.yml
conda activate epig4nn
```

## G4 formation prediction with pre-trained models

This repository contains models trained on various epigenetic marks in A549 cells (see [preprint](https://doi.org/10.1101/2023.03.28.534555) for more details) that can be used for predictions.

To make predictions, python script `predict.py` needs a `.fa` or fasta-style `.txt` file with 1000nt sequences centered at G4, and the respective normalized (0.0-1.0 range) epigenetic arrays in a `.csv` file.
The epigenetic data table should be structured as follows: (n_rows, n_columns) = (n_samples, 1000). The first column should correspond to the sample IDs from the `.fa` file. See `pqs_chr22.txt` and `h3k4me3_293T_chr22.csv` in the `data` folder for an example.
If the prediction is made using a sequence-only model (no epigenetic data), only the sequence file is required.

`predict.py` has the following arguments:
| Argument | Default value | Description
| --- | --- | --- |
| `--seqs` | None | (required) input `.fa` file with 1000nt G4 sequences |
| `--epi` | None | (required) input `.csv` file with 1000nt epigenetic arrays |
| `--model`  | None | (required) model to use: seq, h3k4me3, atac, h3k27ac, h3k4me1, h3k36me3, h3k9me3 |
| `--output` | None | (required) output filename prefix; will be saved to output folder |

Examples of usage:

```bash
#H3K4me3-based model prediction
python predict.py --seqs data/pqs_chr22.txt --epi data/h3k4me3_293T_chr22.csv --model h3k4me3 --output out_h3k4me3

#sequence-based model prediction
python predict.py --seqs data/pqs_chr22.txt --model seq --output out_seq
```

The output will be saved to the `output` folder with the specified prefix.

## Data preparation pipeline for new model training

`data_prep.sh` contains the full data preparation pipeline for new model training. Prerequisites: human genome file, `.bedgraph` with experimental G4 scores, `.bedgraph` with epigenetic scores. A toy dataset for pipeline execution is stored in the `data/tiny` folder. Usage:

1. Download and unpack the human genome from the UCSC website into the `data` folder:

```bash
cd data/
wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz' -O hg19.fa.gz
gunzip hg19.fa.gz
cd ..
```

epiG4NN data preparation pipeline requires `bedtools`. To install, run:

```bash
#ubuntu
apt-get install bedtools

#mac
brew tap homebrew/science
brew install bedtools
```

2. Run the data preparation script with the toy G4 and epigenetic data in `data/tiny`. PQS sequences found in hg19 are used by default.
If you are preprocessing your own data, edit lines 9-10 in `data_prep.sh` to point to your files.

```bash
bash data_prep.sh
```

This script will save the preprocessed `.npy` inputs into `data/inputs_numpy/train_1000nt`.

## Training of new models

After generating the inputs with the `data_prep.sh` script, one can train either a bare sequence model (`train_seq.py`) or a model with epigenetic inputs (`train_epi.py`).

`train_seq.py` and `train_epi.py` have the following arguments:
| Argument | Default value | Description
| --- | --- | --- |
| `--model_name` | None | (required) model name prefix |
| `--inputs` | None | (required) folder with inputs split by chromosome |
| `--kernel_size1` | 11 | kernel size 1 for convolution |
| `--kernel_size2` | 11 | kernel size 2 for convolution |
| `--no_stacks`  | 2 | number of RB stacks (convolutional towers) |
| `--dilation1` | 4 | dilation rate 1 |
| `--dilation2` | 4 | dilation rate 2 |
| `--batch` | 32 | batch size for training |
| `--const_lr` | False | use constant learning rate |

Examples of training:
	
```bash
#H3K4me3-based model training
python train_epi.py --inputs data/inputs_numpy/train_1000nt --model_name tiny_293T_h3k4me3

#sequence-based model training
python train_seq.py --inputs data/inputs_numpy/train_1000nt --model_name tiny_293T_seq
```

## Utils
`PQS_search.py` file contains a python re implementation for PQS search in the human genome

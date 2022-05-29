#! /bin/sh


stage=0
stop_stage=2

. ./parse_options.sh || exit 1;

# 0. Sorting the data files

if [ ${stage} -le 0 ] && [ ${stop_stage} -ge 0 ]; then


    sort \
        -k1,1 \
        -k2,2n \
        K562_H3K27ac.bgr \
        > k27ac_sorted.bgr

fi

# 1. Take intersection of G4 scores, methylation and nucleosome depletion files with filtered PQS file

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then


    bedtools intersect \
        -a PQS_all_merged_padded_250_singleID.bed.sorted \
        -b k27ac_sorted.bgr \
        -wa -wb \
        -sorted \
        > k27ac_intersect.bed
fi

# 2. Merge the records in scores, methylation and nucleosome depletion files

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then

    #cd ../code

    python3 merge_data.py \
        --method meth \
        --read k27ac_intersect.bed \
        --write k27ac_final.bed 
fi

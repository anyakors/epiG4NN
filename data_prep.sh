#! /bin/sh
# contributed by Jasraj Singh @ NTU

stage=0
stop_stage=3

. ./parse_options.sh || exit 1;

epigen_data=data/tiny/ENCFF315TAU_chr22.bedgraph
g4_scores_file=data/tiny/GSE133379_chr22.bedgraph
hg19=data/hg19.fa

# 0. Sorting the data files

if [ ${stage} -le 0 ] && [ ${stop_stage} -ge 0 ]; then

    sort \
        -k1,1 \
        -k2,2n \
        $g4_scores_file \
        > data/g4_scores_sorted.bgr

    sort \
        -k1,1 \
        -k2,2n \
        $epigen_data \
        > data/epi_sorted.bgr

fi

# 1. Take intersection of G4 scores, methylation and nucleosome depletion files with filtered PQS file

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then

    bedtools intersect \
        -a data/PQS_unpadded.bed \
        -b data/g4_scores_sorted.bgr \
        -wa -wb \
        -sorted \
        > data/g4_scores_intersect.bed

    bedtools intersect \
        -a data/PQS_padded.bed \
        -b data/epi_sorted.bgr \
        -wa -wb \
        -sorted \
        > data/epi_intersect.bed
fi

# 2. Merge the records in scores, methylation and nucleosome depletion files

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then

    python3 data/merge_data.py \
        --method scores \
        --read g4_scores_intersect.bed \
        --write g4_scores_final.bed 

    python3 data/merge_data.py \
        --method meth \
        --read epi_intersect.bed \
        --write epi_final.bed 
fi


# 8. Prepare input from the generated files

if [ ${stage} -le 3 ] && [ ${stop_stage} -ge 3 ]; then

    mkdir data/epi_bychr

    python data/normalize_bed.py

    for chr in `cut -f 1 data/epi_final_norm.bed | sort | uniq`; do
        grep -w $chr data/epi_final_norm.bed > data/epi_bychr/$chr.bed
    done

    python data/save_inputs.py \
        --hg19 $hg19 \

    rm epi_final* epi_intersect.bed epi_sorted.bgr g4_scores_final.bed g4_scores_intersect.bed g4_scores_sorted.bgr

fi

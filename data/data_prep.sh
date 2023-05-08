#! /bin/sh
# contributed by Jasraj Singh @ NTU

stage=0
stop_stage=3

. ./parse_options.sh || exit 1;

epigen_data=GSE133379_293T-G4P-hg19-subtract-rep1.bgr
g4_scores_file=GSE133379_293T-G4P-hg19-subtract-rep1.bgr

# 0. Sorting the data files

if [ ${stage} -le 0 ] && [ ${stop_stage} -ge 0 ]; then

    sort \
        -k1,1 \
        -k2,2n \
        $g4_scores_file \
        > g4_scores_sorted.bgr

    sort \
        -k1,1 \
        -k2,2n \
        $epigen_data \
        > epi_sorted.bgr

fi

# 1. Take intersection of G4 scores, methylation and nucleosome depletion files with filtered PQS file

if [ ${stage} -le 1 ] && [ ${stop_stage} -ge 1 ]; then

    bedtools intersect \
        -a PQS_unpadded.bed \
        -b g4_scores_sorted.bgr \
        -wa -wb \
        -sorted \
        > g4_scores_intersect.bed

    bedtools intersect \
        -a PQS_padded.bed \
        -b epi_sorted.bgr \
        -wa -wb \
        -sorted \
        > epi_intersect.bed
fi

# 2. Merge the records in scores, methylation and nucleosome depletion files

if [ ${stage} -le 2 ] && [ ${stop_stage} -ge 2 ]; then

    python3 merge_data.py \
        --method scores \
        --read g4_scores_intersect.bed \
        --write g4_scores_final.bed 

    python3 merge_data.py \
        --method meth \
        --read epi_intersect.bed \
        --write epi_final.bed 
fi


# 8. Prepare input from the generated files

if [ ${stage} -le 3 ] && [ ${stop_stage} -ge 3 ]; then

    mkdir epi_bychr

    python normalize_bed.py

    for chr in `cut -f 1 epi_final_norm.bed | sort | uniq`; do
        grep -w $chr epi_final_norm.bed > epi_bychr/$chr.bed
    done

    python save_inputs.py \
        --hg19 hg19.fa \

fi

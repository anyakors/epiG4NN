for chr in `cut -f 1 k27ac_final_norm.bed | sort | uniq`; do
    grep -w $chr k27ac_final_norm.bed > k27_final_byChr/$chr.bed

done

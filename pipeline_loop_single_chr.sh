#!/bin/bash

#### only for single chromosome !!!

dir_code=$1 
file_pred_site=$2
file_RAWexpected=$3
file_RAWobserved_dict=$4
file_domain=$5
dir_work=$6
chr=$7
mkdir -p $dir_work/candidate_by_chr/
file_random_loop=$dir_work/tmp_random_loop_$chr
file_loop_withindomain=$dir_work/tmp_loop_withindomain_$chr
file_loop_withindomain_withreads=$dir_work/candidate_by_chr/candidate.bed_$chr


grep -w $chr $file_pred_site | sort -k2n,2 | uniq > $dir_work/tmp_pred_site_$chr
python $dir_code/pred2randompair.py $dir_work/tmp_pred_site_$chr $file_random_loop
awk '{print $1"\t"$2"\t"$6"\t"$0}' $file_random_loop | intersectBed -a - -b $file_domain -f 0.99 -wa -wb | cut -f 4-9 > $file_loop_withindomain
python $dir_code/extract_hic_reads_for_loops2.py $file_RAWexpected $file_RAWobserved_dict $file_loop_withindomain $file_loop_withindomain_withreads 



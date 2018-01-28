#!/bin/bash

#### only for single chromosome !!!

dir_code=$(cd `dirname $0`; pwd)
file_pred_site=$1
file_fft=$2
dir_hic_reads=$3
#file_RAWobserved_dict_norm=$3
file_domain=$4
dir_work=$5
chr=$6
mkdir -p $dir_work/candidate_by_chr/
file_random_loop=$dir_work/tmp_random_loop_$chr
file_loop_withindomain=$dir_work/tmp_loop_withindomain_$chr
file_loop_withindomain_withreads=$dir_work/candidate_by_chr/candidate.bed_${chr}.withreads


#bash ~/project1/CISD-master/hicnorm2.sh $dir_hic_reads $chr



file_RAWobserved_dict_norm=$dir_hic_reads/${chr}_5kb.RAWobserved.dict.normalized
grep -w $chr $file_pred_site | sort -k2n,2 | uniq > $dir_work/tmp_pred_site_$chr
python $dir_code/pred2randompair.py $dir_work/tmp_pred_site_$chr $file_random_loop


awk '{print $1"\t"$2"\t"$6"\t"$0}' $file_random_loop | intersectBed -a - -b $file_domain -f 0.99 -wa -wb | cut -f 4-9 > $file_loop_withindomain
#python $dir_code/extract_hic_reads_for_loops2.py $file_RAWexpected $file_RAWobserved_dict $file_loop_withindomain $file_loop_withindomain_withreads 
python $dir_code/extract_hic_reads_for_loops5.py  $file_RAWobserved_dict_norm $file_loop_withindomain $file_loop_withindomain_withreads 


cut -f 1-3 $file_loop_withindomain > $dir_work/tmp1 ; 
cut -f 4-6 $file_loop_withindomain > $dir_work/tmp2; 
Rscript $dir_code/site_by_chrTofft1to50_by_chr_tmp.R $file_fft $dir_work/tmp1  $dir_work/tmp1.fft1to50; 
Rscript $dir_code/site_by_chrTofft1to50_by_chr_tmp.R $file_fft $dir_work/tmp2  $dir_work/tmp2.fft1to50; 
paste $file_loop_withindomain_withreads $dir_work/tmp1.fft1to50 $dir_work/tmp2.fft1to50 | awk '{print $0"\t"int(($6+$5)/2)-int(($3+$2)/2)}' > ${file_loop_withindomain_withreads}.withfft.withdistance; 
rm $dir_work/tmp*





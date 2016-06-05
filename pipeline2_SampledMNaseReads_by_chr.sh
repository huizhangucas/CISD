#!/bin/bash

dir_code=$1
hspeaks=$2
dir_out=$3
dir_fft=$4
model1_SVM=$5
chr=$6

mkdir -p $dir_out
hssites=$dir_out/candidate_${chr}.bed
perl $dir_code/hspeak2points.pl $hspeaks $hssites

######### candidate_to_fft1-50 #####
candidate_sites=$hssites

mkdir -p $dir_out/candidate_by_chr
grep -w $chr $candidate_sites > $dir_out/candidate_by_chr/candidate_${chr}.bed
Rscript $dir_code/site_by_chrTofft1to50_by_chr.R $dir_code $dir_fft $dir_out/candidate_by_chr candidate $chr
cp -f $dir_out/candidate_by_chr/candidate_${chr}.bed.f1to50 $dir_out/

#intersectBed -a $dir_out/candidate_chrAll.bed.f1to50 -b $chia_cobinding -wa | sort -u | sort -R | head -10000 > $dir_out/p_chrAll.bed.f1to50
#intersectBed -a $dir_out/candidate_chrAll.bed.f1to50 -b $chia_all -wa -v | sort -u | sort -R | head -10000 > $dir_out/n_chrAll.bed.f1to50
#cat $p2_f1to50 > $dir_out/p_chrAll.bed.f1to50
#cat $n2_f1to50 > $dir_out/n_chrAll.bed.f1to50

Rscript $dir_code/training2_pred2_by_chr.R $dir_code $model1_SVM $dir_out $chr

mergeBed -i $dir_out/pred2_${chr} > $dir_out/pred2${chr}.merge

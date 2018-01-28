#!/bin/bash

dir=$1 # the folder which contain chr*_5kb.KRnorm, chr*_5kb.VCnorm, chr*_5kb.RAWobserved, chr*_5kb.RAWexpected
i=$2 # chromosome number, for example, 1 for chr1, 2 for chr2, X for chrX

dir_code=$(cd `dirname $0`; pwd)

#for  i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do echo $i; 


#for i in 21; do
filekr=$dir/chr${i}_5kb.KRnorm; 
filevc=$dir/chr${i}_5kb.VCnorm; 
file_obs=$dir/chr${i}_5kb.RAWobserved
file_exp=$dir/chr${i}_5kb.RAWexpected
paste $filekr $filevc | awk '{if($1 !="NaN"){print $1} else {print $2}}' > ${filekr}.tmp; 
python $dir_code/hicnorm.py ${filekr}.tmp $file_exp $file_obs ${file_obs}.dict.normalized 
rm ${filekr}.tmp

#done

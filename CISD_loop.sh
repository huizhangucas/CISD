#!/bin/bash
# This scripts is the main pipeline of CISD_loop. It uses the CISD sites and Hi-C reads as input, according to the selected chromosomes, gives the predicted loops. Before use run this script, you shold run CISD and get the CISD sites in advance. The first parameter of CISD_loop is the output directory of CISD, which is the third parameter of CISD . The second paremeter of CISD_loop is the directory of Hi-C contact matrix. The third parameter of CISD_loop is the chromosomes you want to choose and the fourth parameter of CISD_loop is the output directory of CISD_loop. The usage of this scripts is simple as the following command:
# Usage: bash CISD_loop.sh Input1 Input2 Input3 Output
# Input1: The CISD output directory, which is the same as the third parameter of CISD.
# Input2: The directory containing Hi-C contact matrix and expected reads, which are raw observed reads and raw expected reads. The Hi-C contact matrix must be the same format as the format given in http://dx.doi.org/10.1016/j.cell.2014.11.021. The contact matrix must be at 5kb resolution. The filename of must be the same as the name given in http://dx.doi.org/10.1016/j.cell.2014.11.021, like chr1_5kb.RAWobserved and chr1_5kb.RAWexpected. For more information, please refer to http://dx.doi.org/10.1016/j.cell.2014.11.021.

# Input3: The chromosomes that you want to choose. Different chromosomes should be seperated by "," and it is strongly recommended to do the prediction on all chromosomes. For example, if you want to do the prediction on chromosome 1 and chromosome 2, you may set this parameter as 1,2. If you want to do the prediction on all chromosomes of human, you may set this parameter as 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X.
# Output: The output directory of CISD, the predicted loops will be saved as CISD_loop.txt in this directory.

dir_cisd_site=$1    
dir_hic_reads=$2;   
chr_list=$3;        
dir_out2=$4;        
dir_code=$(cd `dirname $0`; pwd);
mkdir -p $dir_out2/tmp2/;
echo "running CISD_loop......"
#for i in ${chr_list//,/ } ; do awk '{print $1"_"$2"\t"$3}' ${dir_hic_reads}/chr${i}_5kb.RAWobserved > $dir_out2/tmp2/chr${i}_5kb.RAWobserved.dict; echo chr$i; bash ${dir_code}/pipeline_loop_single_chr.sh ${dir_code} $dir_cisd_site/high_score_peaks/hspeaks0.50_pred2/pred2chr${i}.merge  ${dir_hic_reads}/chr${i}_5kb.RAWexpected  $dir_out2/tmp2/chr${i}_5kb.RAWobserved.dict  ${dir_code}/data/domain.bed   ${dir_out2} chr${i};done;


for i in ${chr_list//,/ } ; do grep -w chr$i $dir_cisd_site/CISD_site.txt > $dir_cisd_site/tmp.txt ; awk '{print $1"_"$2"\t"$3}' ${dir_hic_reads}/chr${i}_5kb.RAWobserved > $dir_out2/tmp2/chr${i}_5kb.RAWobserved.dict; echo chr$i; bash ${dir_code}/pipeline_loop_single_chr.sh ${dir_code} $dir_cisd_site/tmp.txt  ${dir_hic_reads}/chr${i}_5kb.RAWexpected  $dir_out2/tmp2/chr${i}_5kb.RAWobserved.dict  ${dir_code}/data/domain.bed   ${dir_out2} chr${i};done;

cat $dir_out2/candidate_by_chr/candidate.bed_chr* | sort -k1.4n,1 -k2n,2 | uniq > $dir_out2/candidate.bed_chrAll;

Rscript ${dir_code}/pred_loop_svm.R  ${dir_code}/data/model2_SVM  $dir_out2/candidate.bed_chrAll  $dir_out2/CISD_loop.txt;

echo "Your CISD loops are saved to $dir_out2/CISD_loop.txt!";

rm -rf $dir_out2/tmp*;


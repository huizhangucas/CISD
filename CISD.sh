#!/bin/bash
# This scripts is the main pipeline of CISD. It uses the MNase-seq signal as input, according to the selected chromosomes, automatically normalize the data, calculate FFT profiles, call high score peaks and finally gives the CISD sites. Before use run this script, you shold run iNPS to denoise the MNase-seq signal. The -o parameter of iNPS is the first parameter of CISD. The second paremeter of CISD is the chromosomes you want to choose and the third parameter of CISD is the output directory of CISD. The usage of this scripts is simple as the following command:
# Usage: bash CISD.sh Input1 Input2 Output
# Input1: The iNPS output directory with prefix, which is the same as the -o parameter of iNPS. CISD uses the .like_wig format data, which is in the iNPS output directory. For more information, please refer to http://dx.doi.org/10.1038/ncomms5909.
# Input2: The chromosomes that you want to choose. Different chromosomes should be seperated by "," and it is strongly recommended to do the prediction on all chromosomes. For example, if you want to do the prediction on chromosome 1 and chromosome 2, you may set this parameter as 1,2. If you want to do the prediction on all chromosomes of human, you may set this parameter as 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X.
# Output: The output directory of CISD, the predicted CISD sites will be saved as CISD_site.txt in this directory.

dir_inps_prefix=$1; 
chr_list=$2;        
dir_out=$3;         
dir_code=$(cd `dirname $0`; pwd);
mkdir -p $dir_out/wig;
mkdir -p $dir_out/tmp;

echo -e "calculating the standard deviation......";
for i in ${chr_list//,/ } ; do sed '1,3d' ${dir_inps_prefix}_chr${i}.like_wig | awk '{print $4}' > $dir_out/wig/chr${i}_4;echo chr$i; done;
cat $dir_out/wig/chr*_4 > $dir_out/tmp/chrAll_4;

sd=$(Rscript ${dir_code}/wig2Sd.R $dir_out/tmp/chrAll_4);
echo -e "standard deviation:${sd}\nnormalizing......";
for i in ${chr_list//,/ } ; do sed '1,3d' ${dir_inps_prefix}_chr${i}.like_wig  |awk '{printf "%.4f\n", $4/"'"$sd"'"}' > $dir_out/wig/chr${i}_4.normalized;echo chr$i; done;

for i in ${chr_list//,/ } ; do Rscript ${dir_code}/calculate_fft_by_chr.R ${dir_code} $dir_out/wig/ $dir_out/wig/  chr${i}; done;

for i in ${chr_list//,/ } ; do Rscript ${dir_code}/calculate_logisticRegression_by_chr.R ${dir_code} $dir_out/wig/ ${dir_code}/data/model1_LRM $dir_out/wig/ chr${i}; done;

mkdir -p $dir_out/high_score_peaks/hspeaks0.50_by_chr/;

for i in ${chr_list//,/ } ; do Rscript ${dir_code}/get_hspeaks_by_chr.R ${dir_code} 0.50 $dir_out/wig/ $dir_out/high_score_peaks/hspeaks0.50_by_chr/ chr${i}; done;

cat $dir_out/high_score_peaks/hspeaks0.50_by_chr/hspeaks0.50chr*.bed > $dir_out/high_score_peaks/hspeaks0.50chrAll.bed

for i in ${chr_list//,/ } ; do bash ${dir_code}/pipeline2_SampledMNaseReads_by_chr.sh ${dir_code} $dir_out/high_score_peaks/hspeaks0.50_by_chr/hspeaks0.50chr${i}.bed $dir_out/high_score_peaks/hspeaks0.50_pred2 $dir_out/wig/ ${dir_code}/data/model1_SVM chr${i};done;

cat $dir_out/high_score_peaks/hspeaks0.50_pred2/pred2chr*.merge > $dir_out/CISD_site.txt;

echo "Your CISD sites are saved to $dir_out/CISD_site.txt!";

rm -rf $dir_out/tmp*;


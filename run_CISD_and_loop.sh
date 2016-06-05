#!/bin/bash

dir_inps_prefix=$1; # Input, which is the same as the -o parameter of iNPS. 
                    # The files that used in this folder must be .like_wig format.
                    # For more information, please refer to http://dx.doi.org/10.1038/ncomms5909

dir_hic_reads=$2;   # Input, which is the folder containing the Hi-C contact matrix and expected reads, which are raw observed reads and raw expected reads.
                    # The Hi-C contact matrix must be the same format as the format given in http://dx.doi.org/10.1016/j.cell.2014.11.021; 
                    # The contact matrix must be at 5kb resolution; 
                    # The filename of must be the same as the name given in http://dx.doi.org/10.1016/j.cell.2014.11.021, like chr1_5kb.RAWobserved and chr1_5kb.RAWexpected.
                    # For more information, please refer to http://dx.doi.org/10.1016/j.cell.2014.11.021

chr_list=$3;        # Input, which are the chromosomes that you want to do prediction.
                    # The format must be like XX,XX,XX,XX XX refers to different chromosomes. For example, 1,2,X means doing prediction on chromosome1, chromosome2 and chromosomeX. 
                    # If you want to do the prediction on whole genome, please use 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X
                    # It is strongly to do the prediction on the whole genome rather than a few chromosomes in one time.
                    # For each choosed chromosome, you must have the corresponding .like_wig file, the .RAWobserved file and the .RAWexpected in the directories mentioned above.

dir_out=$4;         # CISD and CISD_loop output directory.
                    # The predicted CISD sites and loops will be saved as CISD_site.txt and CISD_loop.txt in this directory.




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

dir_cisd_site=$dir_out
dir_out2=$dir_out

mkdir -p $dir_out2/tmp2/;

echo "running CISD_loop......"
for i in ${chr_list//,/ } ; do awk '{print $1"_"$2"\t"$3}' ${dir_hic_reads}/chr${i}_5kb.RAWobserved > $dir_out2/tmp2/chr${i}_5kb.RAWobserved.dict; echo chr$i; bash ${dir_code}/pipeline_loop_single_chr.sh ${dir_code} $dir_cisd_site/high_score_peaks/hspeaks0.50_pred2/pred2chr${i}.merge  ${dir_hic_reads}/chr${i}_5kb.RAWexpected  $dir_out2/tmp2/chr${i}_5kb.RAWobserved.dict  ${dir_code}/data/domain.bed   ${dir_out2} chr${i};done;

cat $dir_out2/candidate_by_chr/candidate.bed_chr* | sort -k1.4n,1 -k2n,2 | uniq > $dir_out2/candidate.bed_chrAll;

Rscript ${dir_code}/pred_loop_svm.R  ${dir_code}/data/model2_SVM  $dir_out2/candidate.bed_chrAll  $dir_out2/CISD_loop.txt;

echo "Your CISD loops are saved to $dir_out2/CISD_loop.txt!";
#rm -rf $dir_out/tmp*;
#rm -rf $dir_out2/tmp*;


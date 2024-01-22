#!/bin/bash
root_dir=`pwd`
samplist=$1 
f=$2
b=$3
date=$4
shell=$root_dir/scripts/other/sumpos/shell2/

perlcombine=/date/guoyuqing/pip/MCTA/perl/sucombine.pl

out1_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/sumT_C_reads
out2_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/sumT_C_mepm
out3_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/sumP_C_reads
out4_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/sumP_C_umi
mkdir -p $out1_dir $out2_dir $out3_dir $out4_dir 


out1=$out1_dir/sumT_${date}_C_reads.txt
out2=$out2_dir/sumT_${date}_C_mepm.txt
out3=$out3_dir/sumP_${date}_C_reads.txt
out4=$out4_dir/sumP_${date}_C_umepm.txt

sumtcombine=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/sumTcombine
sumpcombine=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/sumPcombine


echo "
export PERL5LIB=/home/software/perl5/lib/perl5::/home/software/perl5/lib/perl5:/home/software/perl5/lib/perl5/x86_64-linux-thread-multi

perl $perlcombine $sumtcombine $samplist 7 8 $out1
perl $perlcombine $sumtcombine $samplist 7 9 $out2
perl $perlcombine $sumpcombine $samplist 7 8 $out3
perl $perlcombine $sumpcombine $samplist 7 9 $out4
" > $shell/sumpileup_step2_${date}_${f}_${b}.sh

qsub -e $shell -o $shell -cwd -l vf=2g,p=1  $shell/sumpileup_step2_${date}_${f}_${b}.sh


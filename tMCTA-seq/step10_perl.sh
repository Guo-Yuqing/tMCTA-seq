#!/bin/bash
root_dir=`pwd`
perl=/date/guoyuqing/MCTA/pip/perl/combine_new_v2.pl
perlcombine=/date/guoyuqing/MCTA/pip/perl/sucombine.pl
dir=$root_dir/project/analysis/1.normal
tcombine=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/Tcombine
pcombine=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/Pcombine

samplist=$1 
f=$2
b=$3
date=$4
shell=$root_dir/scripts/normal/shell2/

out1_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/T_C_reads
out2_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/T_CT_reads
out3_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/T_C_mepm
out4_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/T_CT_mepm
out5_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/P_C_reads
out6_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/P_CT_reads
out7_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/P_C_umi
out8_dir=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/P_CT_umi

mkdir -p $tcombine $pcombine  $shell $out1_dir $out2_dir $out3_dir $out4_dir $out5_dir $out6_dir $out7_dir $out8_dir

out1=$out1_dir/T_${date}_C_reads.txt
out2=$out2_dir/T_${date}_CT_reads.txt
out3=$out3_dir/T_${date}_C_mepm.txt
out4=$out4_dir/T_${date}_CT_mepm.txt
out5=$out5_dir/P_${date}_C_reads.txt
out6=$out6_dir/P_${date}_CT_reads.txt
out7=$out7_dir/P_${date}_C_umepm.txt
out8=$out8_dir/P_${date}_CT_umepm.txt



echo "
export PERL5LIB=/home/software/perl5/lib/perl5::/home/software/perl5/lib/perl5:/home/software/perl5/lib/perl5/x86_64-linux-thread-multi

perl $perl $samplist $dir $f $b $date

perl $perlcombine $tcombine $samplist 7 8 $out1
perl $perlcombine $tcombine $samplist 7 9 $out2
perl $perlcombine $tcombine $samplist 7 10 $out3
perl $perlcombine $tcombine $samplist 7 11 $out4
perl $perlcombine $pcombine $samplist 7 8 $out5
perl $perlcombine $pcombine $samplist 7 9 $out6
perl $perlcombine $pcombine $samplist 7 10 $out7
perl $perlcombine $pcombine $samplist 7 11 $out8

" > $shell/step10_perlnew_v2_${date}_${f}_${b}.sh
qsub -e $shell -o $shell -cwd -l vf=2g,p=1 $shell/step10_perlnew_v2_${date}_${f}_${b}.sh




#!/bin/bash
root_dir=`pwd`
samdir=$1 
f=$2
b=$3
date=$4
cat $samdir|grep -v "^#"|while read line
do
 echo "$line">./${date}_temp.txt
 i=`awk '{print $1}' ./${date}_temp.txt`
 s=`awk '{print $2}' ./${date}_temp.txt`
 i1=`awk '{print $3}' ./${date}_temp.txt`
 rm ./${date}_temp.txt
 shell=./script/other/sumpos/shell2/ 
mkdir -p $shell
reads_umi_count=$root_dir/project/analysis/1.normal/${f}_${b}trim/ss_umicount/${i1}
mkdir -p $reads_umi_count
 echo "
export PERL5LIB=/home/software/perl5/lib/perl5::/home/software/perl5/lib/perl5:/home/software/perl5/lib/perl5/x86_64-linux-thread-multi 
 perl1=sumpileup.pl
 perl2=$root_dir/run_shell/sumpileup_gyq_20220802.pl
 sumperlmepmcg=sumMepm_cg.pl 
 sumperlumepmcg=sumupMepm_cg.pl
 in4v2=$root_dir/project/analysis/1.normal/${f}_${b}trim/s06.cgcgma_pcr/${s}.5bp.cgcgmat.rmd.in4v2.gz
 in5v2=$root_dir/project/analysis/1.normal/${f}_${b}trim/s06.cgcgma_pcr/${s}.5bp.cgcgmat.umirmd.in5v2.gz
 in5v2_2=$root_dir/project/analysis/1.normal/${f}_${b}trim/s06.cgcgma_pcr/${s}.5bp.cgcgmat.umirmd.in5v2_2.gz
 cgcgcggpos=/data1/guoyuqing/database/human/hg19/CGI/CGI_cgcgcgg_pos.txt
 cpgpos=/data1/guoyuqing/database/human/hg19/CGI/CGIcgcgcggCG.txt
 sumTcgcgcgg1=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/sumTcombine/${s}.txt
 sumPcgcgcgg1=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/sumPcombine/${s}.txt
 reads_C=$reads_umi_count/${s}.reads.methC.txt
 reads=$root_dir/project/analysis/1.normal/cal/${s}.${f}_${b}reads.txt
 temp=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/${s}.temp.txt
 perl \$perl1 \$in4v2 \$cgcgcggpos \$cpgpos|perl \$sumperlmepmcg - \$reads > \$sumTcgcgcgg1
 perl \$perl2 \$in5v2_2 \$cgcgcggpos \$cpgpos \$reads_C >\$temp
 perl \$perl1 \$in5v2 \$cgcgcggpos \$cpgpos|perl \$sumperlumepmcg - \$reads > \$sumPcgcgcgg1
 rm \$temp
" > $shell/${s}_${f}_${b}.step1.sh

 qsub -e $shell -o $shell -cwd -l vf=2g,p=1 -V $shell/${s}_${f}_${b}.step1.sh
done

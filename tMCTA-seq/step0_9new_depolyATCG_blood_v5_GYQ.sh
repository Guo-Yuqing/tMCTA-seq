#!/bin/bash
##########################################################
## Author: Yuqing Guo
## Version: 5.0
## 2022-10-11
## tmctaseq pip
##########################################################
root_dir=`pip`
samdir=$1 #cmp Pn1
f=$2 #6
b=$3 #12
date=$4 #20160221
cat $samdir|grep -v "^#"|while read line
do
 echo "$line">./${date}_temp.txt
 i=`awk '{print $1}' ./${date}_temp.txt`
 s=`awk '{print $2}' ./${date}_temp.txt`
 i1=`awk '{print $3}' ./${date}_temp.txt`
 rm ./${date}_temp.txt
 shell=./script/$s
 mkdir -p $shell
 echo "

 export PERL5LIB=/home/software/perl5/lib/perl5::/home/software/perl5/lib/perl5:/home/software/perl5/lib/perl5/x86_64-linux-thread-multi

 i=$i #cmp98-1 #ss
 s=$s #Pgac1 #Pnlc1
 f=$f #6
 b=$b #12
 i1=$i1 #20160221

 dir=$root_dir/project/analysis/1.normal
 QCdir=$root_dir/project
 QC_clean=$root_dir/project/cleandata/QC_cutadapter_fastp
 dirperl=/date/guoyuqing/MCTA/perl
 database=/data1/guoyuqing/database/human
 shell=$shell

 samtools=/data1/guoyuqing/software/samtools-1.3.1_makeinstall/bin/samtools

 bismark=/datg/guoyuqing/software/miniconda3/envs/meth/bin/bismark
 bowtie=/datg/guoyuqing/software/bowtie2-2.3.5.1-linux-x86_64/
 intersectBed=/data1/guoyuqing/software/bedtools2/bin/intersectBed

 genome=\$database/hg19/bismark_only_hg19
 CGI=\$database/hg19/CGI/CGI-information-copy2.bed.txt # 8 more columns than CGI-information-copy.bed # revise in v2
 CGIpos=\$database/hg19/CGI/posCGI-information-copy2.bed.txt
 CGIneg=\$database/hg19/CGI/negCGI-information-copy2.bed.txt
 dbSNP=\$database/hg19/dbSNP/human_build_137/vcf_by_chromosome/CHB/processed_data
 dnafa=\$database/hg19/DNA_fa
 fa=\$database/hg19/bismark_only_hg19/hg19_24.fa
 dCGI=\$database/hg19/CGI/duplicgi.txt
 dCGI2=\$database/hg19/CGI/duplicgi2.txt
 sq=\$database/hg19/samtools_header/sq.txt
 cgcgcggpos=\$database/hg19/CGI/CGI_cgcgcgg_pos.txt
 cpgpos=\$database/hg19/CGI/CGIcgcgcggCG.txt

 perlQC=/date/guoyuqing/MCTA/run/QCnolane.pl
 perlfill=\$dirperl/read1_2_filter_adapter.pl
 perlxy=\$dirperl/rm_firstx_leny.pl
 perlch3=\$dirperl/ch3deleate.pl
 perlnoch1=\$dirperl/covStrBed.pl
 perlnoch2=\$dirperl/modCytSplt.pl
 perlsort=\$dirperl/sorcharacter.pl
 perlcgcgmat=\$dirperl/extractCGx2.pl
 perlrmv1=\$dirperl/rmduppcrv3.pl
 perlrm1=\$dirperl/umitargetv2_gyq_20220626.pl
 perlrm2=$root_dir/run_shell/umitargetv2_gyq_20220802.pl  ###echo umi_reads.txt
 qcperl=\$dirperl/qc5bp.pl
 perlbampick=\$dirperl/bampick.pl
 perlcal1=\$dirperl/cal_readsv2.pl
 perlcal2=\$dirperl/cal_reads.pl
 perlmepm1=\$dirperl/Mepm_v2.pl #revise in v2
 perlmepm2=\$dirperl/upMepm_v2.pl #revise in v2
 perlpile=\$dirperl/pileup.pl
 perlmepmcg=\$dirperl/Mepm_cg.pl
 perlumepmcg=\$dirperl/upMepm_cg.pl
 perlcombine=\$dirperl/sucombine.pl
 perlchangeflag=/data1/wenlu/perl/cmp/changeflag.pl
 perldepolyATCG=/data1/wenlu/perl/cmp/depolyATCG.pl

 out1=\$dir/\${f}_\${b}trim/s04.align/\$s
 out2=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s
 out3=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr
 out4=\$dir/\${f}_\${b}trim/s07.Tcgibed
 out5=\$dir/\${f}_\${b}trim/s07.Pcgibed
 out6=\$QCdir/cleandata/\$i1/\$i/\$i.1.clean.depolyATCG.fq.gz
 out7=\$QCdir/cleandata/\$i1/\$i/\$i.2.clean.depolyATCG.fq.gz


 raw1=\$QCdir/rawdata/\$i1/\$i/\${i}_*1.f*q.gz
 raw2=\$QCdir/rawdata/\$i1/\$i/\${i}_R2.fastq.gz
 R1_trim=\$QCdir/R1_trim_rawdata/\$i1/\$i/\${i}_0trim1.fq.gz
 c1=\$QC_clean/\$i1/\$i/\${i}_clean2.1.fq.gz
 c2=\$QC_clean/\$i1/\$i/\${i}_clean2.2.fq.gz
 dir1=\$QCdir/cleandata/sample/\$s.1.clean.fq.gz
 dir2=\$QCdir/cleandata/sample/\$s.2.clean.fq.gz
 adap1=\$dir/s01.filter/\$s.1.adapter.fq.gz
 adap2=\$dir/s01.filter/\$s.2.adapter.fq.gz
 fill1=\$dir/\${f}_\${b}trim/s02.rm\${f}_\${b}/\$s.1.filter.fq.gz
 fill2=\$dir/\${f}_\${b}trim/s02.rm\${f}_\${b}/\$s.2.filter.fq.gz
 nonCG2=\$dir/\${f}_\${b}trim/s03.\${f}_\${b}nonCG/\$s.2.nonCG.fq.gz
 outbis=\$dir/\${f}_\${b}trim/s04.align/\$s
 bam1=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.1.filter_bismark.bam
 bam2=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.2.nonCG_bismark.bam
 bam3=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.2.filter_bismark.bam
 sort1=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.1.filter_bismark.sort
 sort2=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.2.nonCG_bismark.sort
 sort3=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.2.filter_bismark.sort
 sbam1=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.1.filter_bismark.sort.bam
 sbam2=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.2.nonCG_bismark.sort.bam
 sbam3=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.2.filter_bismark.sort.bam
 outnoch=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s
 cov1=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/\$s
 cov2=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmat/\$s
 cov3=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmatrmd/\$s
 cgcgmatout=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr

 in7=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5qc.gz
 in4=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.gz
 in10=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.cgcgmat.qc.gz
 in5=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.rmd.gz
 in6=\$dir/\${f}_\${b}trim/s07.Pcgibed/\$s.\${b}bp.CGI.bed
 in9=\$dir/\${f}_\${b}trim/s07.Tcgibed/\$s.\${b}bp.CGI.bed
 in11=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.2.filter_bismark_SE_report.txt
 in12=\$dir/\${f}_\${b}trim/s04.align/\$s/\$s.1.filter_bismark_SE_report.txt
 in13=\$dir/${f}_${b}trim/s07.Pcgibed/$s.${b}bp.qc.dCGI.bed
 in14=\$dir/${f}_${b}trim/s07.Pcgibed/$s.${b}bp.qc.dCGI2.bed
 in15=\$dir/${f}_${b}trim/s07.Pcgibed/$s.${b}bp.rmd.dCGI.bed
 in16=\$dir/${f}_${b}trim/s07.Pcgibed/$s.${b}bp.rmd.dCGI2.bed
 in17=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.exact.gz
 in18=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.rmd.exact.gz
 in19=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.posexact.gz
 in20=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.negexact.gz
 in21=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.rmd.posexact.gz
 in22=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.rmd.negexact.gz
 in4v2=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.rmd.in4v2.gz ##in4
 in5v2=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.umirmd.in5v2.gz ##in5
 in5v2_2=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/\$s.5bp.cgcgmat.umirmd.in5v2_2.gz 
 umi_txt=\$dir/\${f}_\${b}trim/s06.cgcgma_pcr/txt/\$s.umi.txt

 tempbamcgcg=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmat/\$s.bam
 sorttempbamcgcg=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmat/\$s.sort
 stempbamcgcg=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmat/\$s.sort.bam
 outs05cgcgmat=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmat
 tempbamrmd=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmatrmd/\$s.bam
 sorttempbamrmd=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmatrmd/\$s.sort
 stempbamrmd=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmatrmd/\$s.sort.bam
 outs05rmd=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmatrmd
 outs07pcgi=\$dir/\${f}_\${b}trim/s07.Pcgibed
 outs07tcgi=\$dir/\${f}_\${b}trim/s07.Tcgibed
 bed1=\$dir/\${f}_\${b}trim/bed/\$s.\${b}bp.cgcgmat.rm.bed
 bed1in5=\$dir/\${f}_\${b}trim/bed/\$s.\${b}bp.cgcgmat.rm.bed.in5
 point1=\$dir/\${f}_\${b}trim/point/\$s.\${b}bp.cgcgmat.rm.point
 point1in5=\$dir/\${f}_\${b}trim/point/\$s.\${b}bp.cgcgmat.rm.point.in5
 bed2=\$dir/\${f}_\${b}trim/bed/\$s.\${b}bp.bed
 qcbed=\$dir/\${f}_\${b}trim/bed/\$s.\${b}bp.cgcgmat.qc.rm.bed
 qcpoint=\$dir/\${f}_\${b}trim/point/\$s.\${b}bp.cgcgmat.qc.rm.point
 point2=\$dir/\${f}_\${b}trim/point/\$s.\${b}bp.point
 qccgi=\$dir/\${f}_\${b}trim/s07.Pcgibed/\$s.\${b}bp.qc.CGI.bed
 reads=\$dir/cal/\$s.\${f}_\${b}reads.txt
 tmepm=\$dir/\${f}_\${b}trim/s09.mepm/T_\$s.\${b}bp.Mepm.bed
 pmepm=\$dir/\${f}_\${b}trim/s09.mepm/P_\$s.\${b}bp.Mepm.bed
 upmepm=\$dir/\${f}_\${b}trim/s09.mepm/uPnor_\$s.\${b}bp.Mepm.bed
 temp_pos1=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/temp_pos1.txt
 temp_pos2=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/temp_pos2.txt
 temp_neg1=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/temp_neg1.txt
 temp_neg2=\$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/temp_neg2.txt
 Tcgcgcgg1=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/Tcombine/\${s}.txt
 Pcgcgcgg1=$root_dir/project/analysis/7.additional_pipeline/06.normalCGIcgcgcggCG/Pcombine/\${s}.txt

################STEP0_mkdir.sh###################################
 mkdir -p \$QCdir \$QCdir/cleandata/sample \$dir/\${f}_\${b}trim/s02.rm\${f}_\${b} \$dir/\${f}_\${b}trim/s03.\${f}_\${b}nonCG/ \$out1 \$dir/\${f}_\${b}trim/s04.align/\$s/\${s}_splchr \$out2 \$out3 \$out4 \$dir/cal \$dir/\${f}_\${b}trim/bed/ \$dir/\${f}_\${b}trim/s08.Tcgibednorm/ \$out5 \$dir/\${f}_\${b}trim/s08.Pcgibednorm/ \$dir/\${f}_\${b}trim/s09.mepm \$dir/\${f}_\${b}trim/point/ \$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmat \$dir/\${f}_\${b}trim/s05.noch_arrange/\$s/cgcgmatrmd  && \
 echo "step0done"
##############STEP0_mkdir.sh END###############################

################STEP_trimR1_beforeUMI.sh#######################
ln -s  \$raw1 \$R1_trim  && \\
echo "R1_ln_done"
################STEP_trimR1_beforeUMI.sh END####################


##############STEP_QC.sh#######################################
 touch \$shell/\$s.QC.run... && \\
 sh $root_dir/run_shell/cutadapter_blood.sh \$i \$i1

 cp \$c1 \$dir1 && \\
 cp \$c2 \$dir2 && \\
 mv \$shell/\$s.QC.run... \$shell/\$s.QC.ok
###############STEP_QC.sh END#################################

#############STEP1_fill_adapt.sh##############################
 perl \$perlxy \$dir2 \$fill2 \$f 20 \$b && \\
 perl \$perlxy \$dir1 \$fill1 \$b 20 \$f && \\
 echo "\$s.rm_\${f}_\${b}_done"
###################chr3delete##################################
 perl \$perlch3 \$fill2 \$nonCG2 && \\
 echo "\${s}.s03-rm_tandem_nonCG_done"
##################ch3delete END##############################
 \$bismark --bowtie1 --non_directional --unmapped --fastq --phred33-quals --samtools_path /data1/guoyuqing/software/samtools-1.3.1_makeinstall/bin --path_to_bowtie \$bowtie --output_dir \$outbis/ --temp_dir \$outbis/ \$genome  \$fill1 && \\
 \$samtools sort -o \$sbam1 \$bam1 && \\
 \$samtools index \$sbam1 && \\
 \$bismark --bowtie1 --non_directional --unmapped --fastq --phred33-quals --samtools_path /data1/guoyuqing/software/samtools-1.3.1_makeinstall/bin --path_to_bowtie \$bowtie --output_dir \$outbis/ --temp_dir \$outbis/ \$genome  \$nonCG2 && \\
 \$samtools sort -o \$sbam2 \$bam2 && \\
 \$samtools index \$sbam2 && \\

 \$bismark --bowtie1 --non_directional --unmapped --fastq --phred33-quals --samtools_path /data1/guoyuqing/software/samtools-1.3.1_makeinstall/bin --path_to_bowtie \$bowtie --output_dir \$outbis/ --temp_dir \$outbis/ \$genome  \$fill2 && \\
 \$samtools sort -o \$sbam3 \$bam3 && \\
 \$samtools index \$sbam3 && \\

 rm \$bam1 \$bam2 \$bam3 && \\
 echo "s04done"

###########STEP6_cgcgma_pcr.sh######################################
 \$samtools view \$sbam1|perl \$perlchangeflag - > \$cgcgmatout/\$s.5temp1 &&\\
 \$samtools view \$sbam2 > \$cgcgmatout/\$s.5temp2 &&\\
 perl \$perlcgcgmat \$cgcgmatout/\$s.5temp2 \$fa \$in4 &&\\
 perl \$qcperl \$dir1 5 \$in7 &&\\
 perl \$perlbampick \$in7 \$in4 \$in10 &&\\
 perl \$perlrmv1 \$in7 \$in4 5 \$cgcgmatout/\$s.5temp1 \$in5 &&\\
 perl \$perlrm1 \$in7 \$in4 5 \$cgcgmatout/\$s.5temp1 \$cgcgcggpos normal \$in4v2 &&\\
 perl \$perlrm1 \$in7 \$in4 5 \$cgcgmatout/\$s.5temp1 \$cgcgcggpos rmd \$in5v2 &&\\
 perl \$perlrm2 \$in7 \$in4 5 \$cgcgmatout/\$s.5temp1 \$cgcgcggpos rmd \$in5v2_2 \$umi_txt &&\\

 rm \$cgcgmatout/\$s.5temp1 \$cgcgmatout/\$s.5temp2 && echo "\${s}s06done"


#####################7.CGI###########################
 zcat \$in5v2|awk '{OFS=\"\\t\";if(and(\$2,16))print \$3,\$4-1,\$4+length(\$10)-1,\$1,\$5,\"-\";else print \$3,\$4-1,\$4+length(\$10)-1,\$1,\$5,\"+\"}' - > \$bed1 && \\
 awk '{if (\$6==\"+\")print \$1\"\\t\"\$2\"\\t\"\$2+1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}{if (\$6==\"-\")print \$1\"\\t\"\$3-1\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}' \$bed1 > \$point1 && \\
 \$intersectBed -a \$CGI -b \$point1 -c > \$outs07pcgi/\$s.\${b}bp.CGI.bed && \\
 rm \$bed1 \$point1

 zcat \$in4v2|awk '{OFS=\"\\t\";if(and(\$2,16))print \$3,\$4-1,\$4+length(\$10)-1,\$1,\$5,\"-\";else print \$3,\$4-1,\$4+length(\$10)-1,\$1,\$5,\"+\"}' - > \$bed2 && \\
 awk '{if (\$6==\"+\")print \$1\"\\t\"\$2\"\\t\"\$2+1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}{if (\$6==\"-\")print \$1\"\\t\"\$3-1\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}' \$bed2 > \$point2 && \\
 \$intersectBed -a \$CGI -b \$point2 -c > \$outs07tcgi/\$s.\${b}bp.CGI.bed && \\
 rm \$bed2 \$point2 && \\
 echo "\$s.\$b.cgi.done"
##############STEP7_Pcgi_distri.sh STEP7_Tcgi_distri.sh END########

############STEP8_Pcalreads.sh#################################
 zcat \$in10|awk '{OFS=\"\\t\";if(and(\$2,16))print \$3,\$4-1,\$4+length(\$10)-1,\$1,\$5,\"-\";else print \$3,\$4-1,\$4+length(\$10)-1,\$1,\$5,\"+\"}' - > \$qcbed && \\
 awk '{if (\$6==\"+\")print \$1\"\\t\"\$2\"\\t\"\$2+1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}{if (\$6==\"-\")print \$1\"\\t\"\$3-1\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}' \$qcbed > \$qcpoint && \\
 \$intersectBed -a \$dCGI -b \$qcpoint -c > \$in13 && \\
 \$intersectBed -a \$dCGI2 -b \$qcpoint -c > \$in14

 zcat \$in5|awk '{OFS=\"\\t\";if(and(\$2,16))print \$3,\$4-1,\$4+length(\$10)-1,\$1,\$5,\"-\";else print \$3,\$4-1,\$4+length(\$10)-1,\$1,\$5,\"+\"}' - > \$bed1in5 && \\
 awk '{if (\$6==\"+\")print \$1\"\\t\"\$2\"\\t\"\$2+1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}{if (\$6==\"-\")print \$1\"\\t\"\$3-1\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}' \$bed1in5> \$point1in5 && \\
 \$intersectBed -a \$dCGI -b \$point1in5 -c > \$in15 && \\
 \$intersectBed -a \$dCGI2 -b \$point1in5 -c > \$in16 && \\
 rm \$qcbed \$qcpoint \$bed1in5 \$point1in5 && \\
 echo "bed\$s.\$b.done"

 perl \$perlcal1 \$dir1 \$s.clean 0.25 > \$dir/cal/\$s.\${f}_\${b}reads.txt
 perl \$perlcal1 \$in7 \$s.qc6.clean 0.25 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 perl \$perlcal1 \$fill2 \$s.filter 0.25 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 perl \$perlcal1 \$nonCG2 \$s.rmfilter 0.25 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 \$samtools view \$sbam2|perl \$perlcal2 - \$s.unique_mapping 1 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 perl \$perlcal1 \$in4v2 \$s.cgcgmat 1 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 perl \$perlcal1 \$in10 \$s.cgcgmat.qc 1 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 perl \$perlcal1 \$in5v2 \$s.cgcgmat.qc.rmd 1 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 awk -v var=\$s 'BEGIN{sum=0}{sum+=\$12}END{printf \"P-CGI-\"var\"\\t\"sum\"\\t\"}' \$in6 >> \$dir/cal/\$s.\${f}_\${b}reads.txt #revise in v2
 awk -v var=\$s 'BEGIN{sum=0}{sum+=\$12}END{printf \"T-CGI-\"var\"\\t\"sum\"\\t\"}' \$in9 >> \$dir/cal/\$s.\${f}_\${b}reads.txt #revise in v2
 awk -v var=\$s 'BEGIN{sum=0}{sum+=\$5}END{printf \"qc-dCGI-\"var\"\\t\"sum\"\\t\"}' \$in13 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 awk -v var=\$s 'BEGIN{sum=0}{sum+=\$5}END{printf \"qc-dCGI2-\"var\"\\t\"sum\"\\t\"}' \$in14 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 awk -v var=\$s 'BEGIN{sum=0}{sum+=\$5}END{printf \"rmd-dCGI-\"var\"\\t\"sum\"\\t\"}' \$in15 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 awk -v var=\$s 'BEGIN{sum=0}{sum+=\$5}END{printf \"rmd-dCGI2-\"var\"\\t\"sum\"\\t\"}' \$in16 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 tail -n 14 \$in11| awk 'OFS=\"\\t\"{printf \$7\"\\t\"\$9\"\\t\"}' -|awk -v var=\$s 'OFS=\"\\t\"{printf \"filter2_nonCG\"var\"\\t\"1-(\$2+\$3)/(\$5+\$6+\$2+\$3)\"\\t\"}' - >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 tail -n 14 \$in12| awk 'OFS=\"\\t\"{printf \$7\"\\t\"\$9\"\\t\"}' -|awk -v var=\$s 'OFS=\"\\t\"{printf \"filter1_nonCG\"var\"\\t\"1-(\$2+\$3)/(\$5+\$6+\$2+\$3)\"\\t\"}' - >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 perl \$perlcal1 \$raw1 \$s.raw 0.25 >> \$dir/cal/\$s.\${f}_\${b}reads.txt
 awk '{printf \$0\"\\n\"}' \$dir/cal/\$s.\${f}_\${b}reads.txt > \$dir/cal/\$s.\${f}_\${b}readstemp.txt #revise in v2
 mv \$dir/cal/\$s.\${f}_\${b}readstemp.txt \$dir/cal/\$s.\${f}_\${b}reads.txt #revise in v2

 echo "s08done"
############STEP8_Pcalreads.sh END##############################

##################STEP9_Mepm.sh####################################
 awk '{print \$12}' \$reads|perl \$perlmepm1 \$in9 - > \$tmepm
 awk '{print \$14}' \$reads|perl \$perlmepm1 \$in6 - > \$pmepm
 perl \$perlmepm2 \$in6 \$reads > \$upmepm
#################STEP9_Mepm.sh END######################

#####################################################
 perl \$perlpile \$in4v2 \$cgcgcggpos \$cpgpos|perl \$perlmepmcg - \$reads > \$Tcgcgcgg1
 perl \$perlpile \$in5v2 \$cgcgcggpos \$cpgpos|perl \$perlumepmcg - \$reads > \$Pcgcgcgg1
 echo "s09done"
" > $shell/$s.${f}_$b.sh
 qsub -o $shell -e $shell -cwd -l vf=20G,p=4 -V $shell/$s.${f}_$b.sh
done

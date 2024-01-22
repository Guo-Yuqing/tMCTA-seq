
for sam in 20231102blood_tMCTA-seq_pancancer
do
   sh /date/guoyuqing/wenlu/run_shell/step0_9new_depolyATCG_blood_v5_GYQ.sh /date/guoyuqing/wenlu/samplelist/${sam}.txt 6 12 $sam
	#sh /date/guoyuqing/wenlu/run_shell/sumpileup_step1v2.blood.sh /date/guoyuqing/wenlu/samplelist/${sam}.txt 6 12 $sam
#	sh /date/guoyuqing/wenlu/run_shell/counts_umi_Reads.sh  /date/guoyuqing/wenlu/samplelist/${sam}.txt 6 12 $sam
  #     sh /date/guoyuqing/wenlu/run_shell/step10_perl.sh /date/guoyuqing/wenlu/samplelist/${sam}.txt 6 12 $sam
 #      sh /date/guoyuqing/wenlu/run_shell/sumpileup_step2.sh /date/guoyuqing/wenlu/samplelist/${sam}.txt 6 12 $sam
  #     sh /date/guoyuqing/wenlu/run_shell/merge_cal.sh /date/guoyuqing/wenlu/samplelist/${sam}.txt 6 12 $sam
done









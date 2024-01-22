#!/bin/bash
samdir=$1 #cmp Pn1
f=$2 #6
b=$3 #12
date=$4 #20160221
cal_merge_dir=/date/guoyuqing/wenlu/project/analysis/1.normal/cal/merge
mkdir -p $cal_merge_dir

cat $samdir|grep -v "^#"|while read line
do
 echo "$line">./${date}_temp.txt
 i=`awk '{print $1}' ./${date}_temp.txt`
 s=`awk '{print $2}' ./${date}_temp.txt`
 i1=`awk '{print $3}' ./${date}_temp.txt`
 rm ./${date}_temp.txt
 cat /date/guoyuqing/wenlu/project/analysis/1.normal/cal/${s}.${f}_${b}reads.txt >> /date/guoyuqing/wenlu/project/analysis/1.normal/cal/merge/${date}_cal.txt
done

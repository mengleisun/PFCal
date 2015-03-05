#!/bin/bash

dir_name=HEPMC

line=(`grep -n "E\ " $1 | cut -f1 -d:|sed -n "1~1000p"`)
line=("${line[@]}" `wc -l $1`)

[ -d "$dir_name" ] || mkdir $dir_name

for i in `seq 0 4`
do
    j=`echo $i+1|bc`
    temp_name=$dir_name/${1}_${i}
    sed -n "1,3p" $1 > $temp_name 
    sed -n "${line[$i]},`echo ${line[$j]}-1|bc`p" $1 >> $temp_name
    sed -n "${line[5]}p" $1 >> $temp_name
done

#!/bin/bash


while getopts s:b:o:a: option
do
case "${option}"
in

s) FILE1=${OPTARG};;
b) FILE2=${OPTARG};;
o) FILE3=${OPTARG};;
a) FILE4=${OPTARG};;

esac
done

for i in $(cat $FILE1); do
    samtools bedcov --reference /path/to/ref $FILE2 /path/to/data/${i}*.bam > ${i}.depth

    echo ${i} > ${i}.depth2
    echo ${i} > ${i}.depth3

    python3 get_cov.py ${i}.depth >> ${i}.depth2
    python3 get_cov2.py ${i}.depth >> ${i}.depth3

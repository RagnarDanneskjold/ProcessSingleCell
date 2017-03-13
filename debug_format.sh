#!/bin/bash

if [ $# -ne 2 ]
then
    echo "usage: myoutfile run_number"
    exit 1
fi

file=$1
num=$2

cut -f5,6 $file | sort -k1,1 > debug_differences/myout_$num

cd debug_differences

Rscript debug.R myout_$num fcout $num

#!/bin/zsh

for ((i=$1;i<=$2;i++))
    do
        python3 "script/test_bencemark.py" > outfile 2>&1 &
    done
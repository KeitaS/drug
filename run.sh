#!/bin/bash

source ~/.bashrc

s=$(($1 * 10))
for i in {0..9}
do
    e=$(($s + $i))
    qsub -l nodes=1:ppn=32 -d $PWD ./sim.sh $e
done

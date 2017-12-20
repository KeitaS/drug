#!/bin/bash

source ~/.bashrc

for i in {0..100}
do
    qsub -l nodes=1:ppn=32 -d $PWD -j oe -v ARG1=$i ./sim.sh
done

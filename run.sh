#!/bin/bash

source ~/.bashrc
for i in {0..100}
do
    python ribo8.py $i &
done

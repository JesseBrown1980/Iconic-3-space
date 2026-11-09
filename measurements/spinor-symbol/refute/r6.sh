#!/bin/bash
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
grep -E "R6|square_to|M\^2|MM|N\^2|NN" $M/4d/MEASURE-4D.hbp | cut -c1-700
echo ===
grep -E "ADD:R6|ADD.R6|R6" $M/8d/MEASURE-8D.hbp | cut -c1-700
echo ===
grep -E "R6" $M/7d/MEASURE-7D.hbp | cut -c1-500

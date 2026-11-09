#!/bin/bash
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
sed -n '39,40p;57p' $M/7d/MEASURE-7D.precompare.hbp | cut -c1-600
grep -E "id=M8\.(G5|J4|K2|H1|B1|J3)\|" $M/8d/MEASURE-8D.hbp | cut -c1-200

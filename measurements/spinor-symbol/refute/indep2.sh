#!/bin/bash
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
cd $M
stat -c '%y %n' 4d/frozen/* 4d/MEASURE-4D-BEFORE-COMPARE.hbp 4d/measure_4d.before-compare.py 7d/MEASURE-7D.precompare.hbp 7d/measure_7d.precompare.py 8d/frozen/* 8d/MEASURED-ONLY-BEFORE-READ-8D.hbp 8d/measure_8d.before-read.py.txt | sort
echo "=== 8D measured block lines 1 to 5 of pre-read"
head -2 8d/MEASURED-ONLY-BEFORE-READ-8D.hbp | cut -c1-900
echo "=== 4D/7D head rows"
head -1 4d/MEASURE-4D-BEFORE-COMPARE.hbp | cut -c1-1200
head -1 7d/MEASURE-7D.precompare.hbp | cut -c1-1200
echo "=== shared non-trivial lines between scripts (len>=50, stripped)"
for a in 4d/measure_4d.py 7d/measure_7d.py 8d/measure_8d.py; do sed 's/^[ \t]*//' $a | awk 'length($0)>=50' | sort -u > /tmp/refute-$(basename $a).lines; done
comm -12 /tmp/refute-measure_4d.py.lines /tmp/refute-measure_7d.py.lines | wc -l
comm -12 /tmp/refute-measure_4d.py.lines /tmp/refute-measure_8d.py.lines | wc -l
comm -12 /tmp/refute-measure_7d.py.lines /tmp/refute-measure_8d.py.lines | wc -l
comm -12 /tmp/refute-measure_7d.py.lines /tmp/refute-measure_8d.py.lines | head -8
comm -12 /tmp/refute-measure_4d.py.lines /tmp/refute-measure_8d.py.lines | head -8
echo "=== final scripts cross refs"
grep -nE "measure/4d|measure/7d|measure/8d|measure_4d|measure_7d|measure_8d|MEASURE-4D|MEASURE-7D|MEASURE-8D|c28f8af6" 4d/measure_4d.py 7d/measure_7d.py 8d/measure_8d.py | grep -vE "^4d/measure_4d.py:.*(measure_4d|MEASURE-4D)|^7d/measure_7d.py:.*(measure_7d|MEASURE-7D)|^8d/measure_8d.py:.*(measure_8d|MEASURE-8D)" | cut -c1-250

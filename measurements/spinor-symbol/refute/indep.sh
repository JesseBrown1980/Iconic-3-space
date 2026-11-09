#!/bin/bash
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
cd $M
echo "=== sidecars of pre-read snapshots"
( cd 4d && sha256sum -c MEASURE-4D-BEFORE-COMPARE.hbp.sha256 ); ( cd 7d && sha256sum -c MEASURE-7D.precompare.hbp.sha256 ); ( cd 8d && sha256sum -c MEASURED-ONLY-BEFORE-READ-8D.hbp.sha256 )
sha256sum 4d/measure_4d.before-compare.py 7d/measure_7d.precompare.py 8d/measure_8d.before-read.py.txt
echo "=== words in pre-read snapshots that name PRED/ADD content"
for f in 4d/MEASURE-4D-BEFORE-COMPARE.hbp 7d/MEASURE-7D.precompare.hbp 8d/MEASURED-ONLY-BEFORE-READ-8D.hbp; do
 echo "--- $f rows=$(wc -l < $f)"
 grep -noiE "addendum[^|]{0,40}|\bADD[ :.][^|]{0,30}|PRED[ :.][^|]{0,30}|predict[^|]{0,40}|R[1-6][_ ][^|]{0,30}|DIMS[^|]{0,20}|GPT[^|]{0,30}" $f | head -20
done
echo "=== pre-read scripts: opens and frozen references"
for f in 4d/measure_4d.before-compare.py 7d/measure_7d.precompare.py 8d/measure_8d.before-read.py.txt; do
 echo "--- $f"
 grep -nE "open\(|argv|frozen|PREDICTION|ADDENDUM|addendum|/4d|/7d|/8d|measure_[478]d|c28f8af6" $f | head -30
done

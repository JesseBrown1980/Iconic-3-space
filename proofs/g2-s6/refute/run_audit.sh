#!/bin/bash
B=/mnt/c/tmp/opus55/g2-s6-20261109
mkdir -p $B/refute/audit/m1 $B/refute/audit/m2 $B/refute/audit/m3
cp $B/m1/m1_octonion.py $B/refute/audit/m1/
cp $B/m2/m2_roots.py $B/refute/audit/m2/
cp $B/m3/m3_towers.py $B/refute/audit/m3/
cd $B/refute
time python3 -I -B $B/refute/audit_float.py > $B/refute/audit_rows.txt 2> $B/refute/audit_err.txt
echo "AUDIT_EXIT=$?"
cat $B/refute/audit_rows.txt
tail -5 $B/refute/audit_err.txt
sha256sum $B/refute/audit/m1/M1-OCTONION.hbp $B/refute/audit/m2/M2-ROOTS.hbp $B/refute/audit/m3/M3-TOWERS.hbp

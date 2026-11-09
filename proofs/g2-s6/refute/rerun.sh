#!/bin/bash
B=/mnt/c/tmp/opus55/g2-s6-20261109
R=$B/refute/rerun
cp $B/m1/m1_octonion.py $R/m1/
cp $B/m2/m2_roots.py $R/m2/
cp $B/m3/m3_towers.py $R/m3/
sha256sum $B/m1/m1_octonion.py $B/m2/m2_roots.py $B/m3/m3_towers.py $R/m1/m1_octonion.py $R/m2/m2_roots.py $R/m3/m3_towers.py
python3 --version
for j in m1:m1_octonion.py m2:m2_roots.py m3:m3_towers.py; do
  d=${j%%:*}; s=${j##*:}
  cd $R/$d
  python3 -I -B $s > out.txt 2> err.txt
  echo "JOB=$d EXIT=$? STDERR_BYTES=$(wc -c < err.txt)"
done
sha256sum $R/m1/M1-OCTONION.hbp $B/m1/M1-OCTONION.hbp $R/m2/M2-ROOTS.hbp $B/m2/M2-ROOTS.hbp $R/m3/M3-TOWERS.hbp $B/m3/M3-TOWERS.hbp
cmp $R/m1/M1-OCTONION.hbp $B/m1/M1-OCTONION.hbp && echo M1_BYTES_EQUAL
cmp $R/m2/M2-ROOTS.hbp $B/m2/M2-ROOTS.hbp && echo M2_BYTES_EQUAL
cmp $R/m3/M3-TOWERS.hbp $B/m3/M3-TOWERS.hbp && echo M3_BYTES_EQUAL
cd $B/m1 && sha256sum -c M1-OCTONION.hbp.sha256
cd $B/m2 && sha256sum -c M2-ROOTS.hbp.sha256
cd $B/m3 && sha256sum -c M3-TOWERS.hbp.sha256
wc -l $R/m*/*.hbp
grep -c -E "FAIL" $R/m1/M1-OCTONION.hbp $R/m2/M2-ROOTS.hbp $R/m3/M3-TOWERS.hbp

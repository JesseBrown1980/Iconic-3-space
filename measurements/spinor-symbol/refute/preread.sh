#!/bin/bash
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
R=/tmp/refute-spinor-preread-20261109
rm -rf $R; mkdir -p $R/4d $R/7d/out $R/8d/out
cp $M/4d/measure_4d.before-compare.py $R/4d/
cp $M/7d/measure_7d.precompare.py $R/7d/
cp $M/8d/measure_8d.before-read.py.txt $R/8d/measure_8d.before-read.py
cd $R/4d && python3 -I -B measure_4d.before-compare.py --measure-only > run.log 2>&1; echo "4d exit=$?"; ls $R/4d
cd $R/7d && python3 -I -B measure_7d.precompare.py $R/7d/out > run.log 2>&1; echo "7d exit=$?"; ls $R/7d/out
cd $R/8d && python3 -I -B measure_8d.before-read.py $R/8d/out > run.log 2>&1; echo "8d exit=$?"; ls $R/8d/out
sha256sum $R/4d/*.hbp $R/7d/out/*.hbp $R/8d/out/*.hbp
echo "--- sealed snapshots"
sha256sum $M/4d/MEASURE-4D-BEFORE-COMPARE.hbp $M/7d/MEASURE-7D.precompare.hbp $M/8d/MEASURED-ONLY-BEFORE-READ-8D.hbp
echo "--- diff 7d"
diff $R/7d/out/MEASURE-7D.hbp $M/7d/MEASURE-7D.precompare.hbp | cut -c1-300 | head
echo "--- diff 8d"
diff $R/8d/out/*.hbp $M/8d/MEASURED-ONLY-BEFORE-READ-8D.hbp | cut -c1-300 | head
echo "--- prefix check: snapshot measured rows inside final file"
python3 -I -B - <<'PY'
M="/mnt/c/tmp/opus55/spinor-symbol-20261109/measure/"
import hashlib
pairs=[("4d/MEASURE-4D-BEFORE-COMPARE.hbp","4d/MEASURE-4D.hbp"),("7d/MEASURE-7D.precompare.hbp","7d/MEASURE-7D.hbp"),("8d/MEASURED-ONLY-BEFORE-READ-8D.hbp","8d/MEASURE-8D.hbp")]
for a,b in pairs:
    A=open(M+a,"rb").read().split(b"\n"); B=open(M+b,"rb").read().split(b"\n")
    n=0
    while n<len(A) and n<len(B) and A[n]==B[n]: n+=1
    Bset=set(B)
    missing=[x[:60] for x in A if x and x not in Bset]
    print(a, "snapshot_rows", len([x for x in A if x]), "common_prefix_rows", n, "snapshot_rows_absent_from_final", len(missing), missing[:4])
PY

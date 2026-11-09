#!/bin/bash
# usage: gate.sh LIST(freeze|transfer) LOG
LIST=$1; LOG=$2
S=/mnt/c/Users/acer/AppData/Local/Temp/claude/C--Users-acer/81118dc5-cc7e-4fd7-a2e3-32a1a175bed9/scratchpad
G=$S/gate-refute-$LIST
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
R=/mnt/c/tmp/opus55/iconic-3-space/repo
rm -rf $G; mkdir -p $G
cd $R
git archive freeze-spinor-symbol-20261109 tools proofs/g2-s6/refute | tar -x -C $G
if [ "$LIST" = "transfer" ]; then
  git archive origin/transfer-theorem-20261109 tools/FLOAT-EVIDENCE.txt proofs/transfer/check proofs/transfer/second-seat | tar -x -C $G
fi
mkdir -p $G/measure/4d $G/measure/7d $G/measure/8d $G/measure/refute
cp $M/4d/measure_4d.py $M/4d/measure_4d.before-compare.py $G/measure/4d/
cp $M/7d/measure_7d.py $M/7d/measure_7d.precompare.py $G/measure/7d/
cp $M/8d/measure_8d.py $M/8d/measure_8d.before-read.py.txt $G/measure/8d/
cp $M/refute/refute_spinor.py $M/refute/refute_spinor.before-compare.py $G/measure/refute/
cd $G
git init -q . && git add -A
{
echo "LIST $LIST"
sha256sum tools/no_float_gate.py tools/FLOAT-EVIDENCE.txt
python3 -I -X utf8 tools/no_float_gate.py --selftest > selftest.out 2>&1; echo "SELFTEST_EXIT $?"; tail -2 selftest.out
python3 -I -X utf8 tools/no_float_gate.py $G > scan.out 2>&1; echo "SCAN_EXIT $?"; tail -5 scan.out
for f in measure/4d/measure_4d.py measure/4d/measure_4d.before-compare.py measure/7d/measure_7d.py measure/7d/measure_7d.precompare.py measure/8d/measure_8d.py measure/8d/measure_8d.before-read.py.txt measure/refute/refute_spinor.py measure/refute/refute_spinor.before-compare.py; do
  python3 -I -X utf8 tools/no_float_gate.py --pins $f 2>&1 | head -2
done
# planted control
printf 'x = 1.5\ny = 7 / 2\n' > plant.py; git add plant.py
python3 -I -X utf8 tools/no_float_gate.py $G > plant.out 2>&1; echo "PLANT_SCAN_EXIT $?"; grep -E "plant.py" plant.out | head -3
git rm -q --cached plant.py; rm plant.py
python3 -I -X utf8 tools/no_float_gate.py $G > scan2.out 2>&1; echo "UNPLANT_SCAN_EXIT $?"
} > $LOG 2>&1
cat $LOG

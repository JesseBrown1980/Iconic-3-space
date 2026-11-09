#!/bin/bash
set -u
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
R=/tmp/refute-spinor-rerun-20261109
rm -rf "$R"; mkdir -p "$R/4d/frozen" "$R/7d/out" "$R/8d/out" "$R/frozen"
cd /mnt/c/tmp/opus55/iconic-3-space/repo
git show freeze-spinor-symbol-20261109:predictions/spinor-symbol/PREDICTION-SPINOR-SYMBOL.hbp > "$R/frozen/PREDICTION-SPINOR-SYMBOL.hbp"
git show freeze-spinor-symbol-20261109:predictions/spinor-symbol/ADDENDUM-SPINOR-SYMBOL.hbp > "$R/frozen/ADDENDUM-SPINOR-SYMBOL.hbp"
sha256sum "$R"/frozen/*
cp "$R"/frozen/* "$R/4d/frozen/"
cp "$M/4d/measure_4d.py" "$R/4d/"
cp "$M/7d/measure_7d.py" "$R/7d/"
cp "$M/8d/measure_8d.py" "$R/8d/"
sha256sum "$R/4d/measure_4d.py" "$R/7d/measure_7d.py" "$R/8d/measure_8d.py"
python3 --version
cd "$R/4d" && ( time python3 -I -B measure_4d.py > run.log 2>&1 ; echo "4d exit=$?" )
cd "$R/7d" && ( time python3 -I -B measure_7d.py "$R/7d/out" "$R/frozen/PREDICTION-SPINOR-SYMBOL.hbp" "$R/frozen/ADDENDUM-SPINOR-SYMBOL.hbp" > run.log 2>&1 ; echo "7d exit=$?" )
cd "$R/8d" && ( time python3 -I -B measure_8d.py "$R/8d/out" "$R/frozen/PREDICTION-SPINOR-SYMBOL.hbp" "$R/frozen/ADDENDUM-SPINOR-SYMBOL.hbp" > run.log 2>&1 ; echo "8d exit=$?" )
echo "--- rerun hashes"
sha256sum "$R/4d/MEASURE-4D.hbp" "$R/7d/out/MEASURE-7D.hbp" "$R/8d/out/MEASURE-8D.hbp"
echo "--- original hashes"
sha256sum "$M/4d/MEASURE-4D.hbp" "$M/7d/MEASURE-7D.hbp" "$M/8d/MEASURE-8D.hbp"
echo "--- cmp"
cmp "$R/4d/MEASURE-4D.hbp" "$M/4d/MEASURE-4D.hbp" && echo 4d_identical
cmp "$R/7d/out/MEASURE-7D.hbp" "$M/7d/MEASURE-7D.hbp" && echo 7d_identical
cmp "$R/8d/out/MEASURE-8D.hbp" "$M/8d/MEASURE-8D.hbp" && echo 8d_identical
cd "$R/4d" && sha256sum -c MEASURE-4D.hbp.sha256
cd "$R/7d/out" && sha256sum -c MEASURE-7D.hbp.sha256
cd "$R/8d/out" && sha256sum -c MEASURE-8D.hbp.sha256
cat "$R/4d/MEASURE-4D.hbp.sha256" "$R/7d/out/MEASURE-7D.hbp.sha256" "$R/8d/out/MEASURE-8D.hbp.sha256"
tail -3 "$R"/*/run.log

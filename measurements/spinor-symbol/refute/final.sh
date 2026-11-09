#!/bin/bash
D=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure/refute
cd $D
echo "--- banned chars in source: slash count, percent count"
grep -c "/" refute_spinor.py; grep -c "%" refute_spinor.py
mkdir -p check && python3 -I -B refute_spinor.py check | head -2
cmp check/REFUTE-SPINOR.measured.hbp REFUTE-SPINOR-BEFORE-COMPARE.hbp && echo MEASURED_BLOCK_UNCHANGED
bash $D/gate.sh freeze $D/GATE-REFUTE-freeze.log > /dev/null 2>&1; tail -12 $D/GATE-REFUTE-freeze.log
bash $D/gate.sh transfer $D/GATE-REFUTE-transfer.log > /dev/null 2>&1; grep -E "EXIT|NO_FLOAT_DONE|refute_spinor.py " $D/GATE-REFUTE-transfer.log
python3 -I -B refute_spinor.py $D full \
  /tmp/refute-spinor-rerun-20261109/frozen/PREDICTION-SPINOR-SYMBOL.hbp \
  /tmp/refute-spinor-rerun-20261109/frozen/ADDENDUM-SPINOR-SYMBOL.hbp \
  $D/REFUTE-SPINOR-BEFORE-COMPARE.hbp \
  /tmp/refute-spinor-rerun-20261109 /tmp/refute-spinor-preread-20261109 \
  /mnt/c/tmp/opus55/spinor-symbol-20261109/measure \
  $D/GATE-REFUTE-freeze.log $D/GATE-REFUTE-transfer.log
echo "exit=$?"
cd $D && sha256sum -c REFUTE-SPINOR.hbp.sha256 && cat REFUTE-SPINOR.hbp.sha256

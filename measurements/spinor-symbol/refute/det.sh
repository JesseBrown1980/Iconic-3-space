#!/bin/bash
D=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure/refute
mkdir -p /tmp/refute-det && cd $D && python3 -I -B refute_spinor.py /tmp/refute-det full \
  /tmp/refute-spinor-rerun-20261109/frozen/PREDICTION-SPINOR-SYMBOL.hbp /tmp/refute-spinor-rerun-20261109/frozen/ADDENDUM-SPINOR-SYMBOL.hbp \
  $D/REFUTE-SPINOR-BEFORE-COMPARE.hbp /tmp/refute-spinor-rerun-20261109 /tmp/refute-spinor-preread-20261109 \
  /mnt/c/tmp/opus55/spinor-symbol-20261109/measure $D/GATE-REFUTE-freeze.log $D/GATE-REFUTE-transfer.log > /dev/null
cmp /tmp/refute-det/REFUTE-SPINOR.hbp $D/REFUTE-SPINOR.hbp && echo DETERMINISTIC
sha256sum $D/GATE-REFUTE-freeze.log $D/GATE-REFUTE-transfer.log $D/REFUTE-SPINOR-BEFORE-COMPARE.hbp $D/refute_spinor.before-compare.py $D/refute_spinor.py
cat $D/REFUTE-SPINOR.hbp.sha256

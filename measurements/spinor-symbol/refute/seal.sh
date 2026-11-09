#!/bin/bash
D=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure/refute
cd $D
cp refute_spinor.py refute_spinor.before-compare.py
cp sealed/REFUTE-SPINOR.measured.hbp REFUTE-SPINOR-BEFORE-COMPARE.hbp
sha256sum REFUTE-SPINOR-BEFORE-COMPARE.hbp | sed 's#  .*#  REFUTE-SPINOR-BEFORE-COMPARE.hbp#' > REFUTE-SPINOR-BEFORE-COMPARE.hbp.sha256
sha256sum -c REFUTE-SPINOR-BEFORE-COMPARE.hbp.sha256
sha256sum refute_spinor.before-compare.py
date -u +%Y-%m-%dT%H:%M:%SZ

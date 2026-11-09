#!/bin/bash
D=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure/refute
cd $D && grep -c "/" refute_spinor.py; python3 -I -B refute_spinor.py sealed && (cd sealed && sha256sum -c REFUTE-SPINOR.measured.hbp.sha256) && cut -c1-420 sealed/REFUTE-SPINOR.measured.hbp

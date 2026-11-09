#!/bin/bash
F=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure/refute/refute_spinor.py
echo "--- slash, percent, double-star, float, f-string scan"
grep -nE "/|%|\*\*|float|f\"|f'" $F | head
cd /mnt/c/tmp/opus55/spinor-symbol-20261109/measure/refute && mkdir -p sealed && timeout 500 python3 -I -B refute_spinor.py sealed; echo "exit=$?"

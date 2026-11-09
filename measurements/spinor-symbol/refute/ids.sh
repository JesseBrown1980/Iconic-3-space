#!/bin/bash
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
grep -E "^MEAS\|id=M4D\.(G1|S1|S7|S4|R1|R2|S8|S9)|^MEAS" $M/4d/MEASURE-4D-BEFORE-COMPARE.hbp | cut -c1-260 | head -45
echo ===
grep -E "^M7" $M/7d/MEASURE-7D.precompare.hbp | grep -E "M0[1-9]|M1[0-9]|R1|R4|assoc|cross|orient|Psi7" | cut -c1-260 | head -30

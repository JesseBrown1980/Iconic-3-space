#!/bin/bash
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
sed -n '78,100p;188,236p' $M/7d/measure_7d.py
echo "=== 8d"
sed -n '64,100p;109,160p' $M/8d/measure_8d.py
echo "=== 4d: how identities compared"
grep -nE "def (ident|poly_identity|check_identity|pm_equal|pm_eq|coeff_check|identity)[a-z_]*\(" $M/4d/measure_4d.py
grep -nE "^def " $M/4d/measure_4d.py | sed -n '1,80p' | tr '\n' ' '

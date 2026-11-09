#!/bin/bash
M=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure
for f in 4d/measure_4d.py 7d/measure_7d.py 8d/measure_8d.py; do
 echo "=== $f"
 grep -nE "^def (p[a-z_]*|poly[a-z_]*|is_zero|check[a-z_]*|holds|cmp[a-z_]*|coef[a-z_]*|mat_eq|eq[a-z_]*|same[a-z_]*|diff[a-z_]*|compare_poly[a-z_]*|count[a-z_]*)\(" $M/$f | head -30
done
echo "=== 8d holds()"
awk '/^def holds\(/,/^def [a-z_]+\(/{print}' $M/8d/measure_8d.py | head -40
echo "=== 7d zero-check helpers"
grep -nE "^def " $M/7d/measure_7d.py | head -60

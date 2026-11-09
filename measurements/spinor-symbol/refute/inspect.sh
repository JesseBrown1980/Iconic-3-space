#!/bin/bash
cd /mnt/c/tmp/opus55/spinor-symbol-20261109/measure
for f in 4d/measure_4d.py 7d/measure_7d.py 8d/measure_8d.py; do
  echo "=== $f"; wc -l "$f"
  grep -nE "sys\.argv|open\(|__file__|os\.path|listdir|chdir|^import|^from|c28f8af6|GPT|PREDICTION-SPINOR|ADDENDUM|def main" "$f" | head -50
done

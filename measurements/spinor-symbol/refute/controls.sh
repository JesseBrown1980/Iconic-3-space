#!/bin/bash
R=/tmp/refute-spinor-rerun-20261109
for f in "$R/4d/MEASURE-4D.hbp" "$R/7d/out/MEASURE-7D.hbp" "$R/8d/out/MEASURE-8D.hbp"; do
  echo "=== $f"
  cut -d'|' -f1 "$f" | sort | uniq -c | tr '\n' ' '; echo
  grep -E "^(CONTROL|CTRL)" "$f" | cut -c1-330
  grep -E "^(MEASURE_SUMMARY|MEASURED_END|MEASURED_BLOCK|COMPARE_SUMMARY|PRE_READ|MEASURED_BLOCK_SEAL)" "$f" | cut -c1-600
done

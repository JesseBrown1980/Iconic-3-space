#!/bin/bash
cd /mnt/c/tmp/opus55/iconic-3-space/repo
git show freeze-spinor-symbol-20261109:predictions/spinor-symbol/PREDICTION-SPINOR-SYMBOL.hbp | grep -E "^(PRED|DIMS|COMPARE|PLANTED|CORRECTION|BOTT)" | cut -c1-700
echo =========== ADDENDUM
git show freeze-spinor-symbol-20261109:predictions/spinor-symbol/ADDENDUM-SPINOR-SYMBOL.hbp

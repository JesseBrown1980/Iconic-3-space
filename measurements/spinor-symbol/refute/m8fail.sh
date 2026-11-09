#!/bin/bash
R=/tmp/refute-spinor-rerun-20261109
f="$R/8d/out/MEASURE-8D.hbp"
grep -E "^M8" "$f" | grep -vE "ok=1|pass=1|holds=1" | cut -c1-900
echo ----
grep -E "^M8" "$f" | grep -nE "ok=0|pass=0|holds=0|FAIL" | cut -c1-900
echo ---- status field values
grep -E "^M8" "$f" | grep -oE "\|(ok|pass|status|result)=[^|]*" | sort | uniq -c

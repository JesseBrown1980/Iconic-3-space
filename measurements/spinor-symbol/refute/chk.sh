#!/bin/bash
F=/mnt/c/tmp/opus55/spinor-symbol-20261109/measure/refute/refute_spinor.py
wc -l $F; grep -c $'\r' $F; awk 'NF==0' $F | wc -l
sed -n '1,4p' $F | od -c | head -5

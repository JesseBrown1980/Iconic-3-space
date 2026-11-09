#!/bin/bash
cd /mnt/c/tmp/opus55/transfer-20261109/proof/check || exit 2
/usr/bin/time -v python3 -I /mnt/c/tmp/opus55/transfer-20261109/proof/check/transfer_check.py /mnt/c/tmp/opus55/transfer-20261109/proof/check/out "$@" 2> /mnt/c/tmp/opus55/transfer-20261109/proof/check/out/time.txt
echo "EXIT=$?"
grep -E "Elapsed|Maximum resident" /mnt/c/tmp/opus55/transfer-20261109/proof/check/out/time.txt

#!/bin/bash
cd /mnt/c/tmp/opus55/iconic-3-space/repo
git show freeze-spinor-symbol-20261109:tools/FLOAT-EVIDENCE.txt
echo ====
git show freeze-spinor-symbol-20261109:tools/no_float_gate.py | sha256sum
git show freeze-spinor-symbol-20261109:tools/FLOAT-EVIDENCE.txt | sha256sum
git show freeze-spinor-symbol-20261109:tools/no_float_gate.py | head -80
echo ....
git show freeze-spinor-symbol-20261109:tools/no_float_gate.py | grep -nE "argv|def main|selftest|git ls-files|subprocess|BANNED|ALLOWED|import" | head -60
git show freeze-spinor-symbol-20261109:tools/no_float_gate.py | wc -l
echo ==== transfer branch list
git show origin/transfer-theorem-20261109:tools/FLOAT-EVIDENCE.txt | grep -vE "^#" | head; git show origin/transfer-theorem-20261109:tools/no_float_gate.py | sha256sum
git show origin/main:tools/no_float_gate.py | sha256sum

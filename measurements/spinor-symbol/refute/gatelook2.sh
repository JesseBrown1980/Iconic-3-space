#!/bin/bash
cd /mnt/c/tmp/opus55/iconic-3-space/repo
git show freeze-spinor-symbol-20261109:tools/no_float_gate.py | sed -n '96,110p'
git show freeze-spinor-symbol-20261109:tools/no_float_gate.py | grep -nE "def main|sys.argv|--pins|--selftest|--file|root =|ROOT" | head -30

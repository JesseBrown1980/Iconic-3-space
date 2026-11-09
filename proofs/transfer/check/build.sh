#!/bin/bash
source /root/.cargo/env 2>/dev/null
cd /mnt/c/tmp/opus55/transfer-20261109/proof/check/enum || exit 2
export CARGO_TARGET_DIR=/root/transfer-enum-target
cargo +1.81 clippy --release -- -D warnings 2>&1 | tail -20
echo "CLIPPY_EXIT=${PIPESTATUS[0]}"
cargo +1.81 build --release 2>&1 | tail -5
echo "BUILD_EXIT=${PIPESTATUS[0]}"
ls -la /root/transfer-enum-target/release/transfer_enum
rustc +1.81 --version

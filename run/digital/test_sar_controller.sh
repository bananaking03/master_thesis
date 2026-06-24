#!/bin/bash
# Run the SAR controller testbench with Cadence Xcelium (`xrun`).
# Mirrors run/digital/test_calibration.sh structure.

# set -euo pipefail

## Resolve repository root directory relative to this script
# SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# : ${WORK_DIR:=$(cd "$SCRIPT_DIR/../.." && pwd)}

xrun \
    -gui \
    -debug \
    -linebug \
    "$WORK_DIR/verilog/tb/sar_adc_conroller_tb.v" \
    "$WORK_DIR/verilog/rtl/sar_controller_hw.v" \

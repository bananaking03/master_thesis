#!/bin/bash
# ============================================================
# Xcelium AMS simulation script for SAR_1024_top
# ============================================================

# --- Paths ---
NETLIST="/users/students/r0884392/simulation/thesis/SAR_un_1024_top_tb/maestro/results/maestro/.tmpADEDir_r0884392/thesis_SAR_un_1024_top_tb_1/thesis_SAR_un_1024_top_tb_schematic_spectre/netlist/input.scs"
VERILOG="../../verilog/rtl/sar_controller_hw.v"
AMS_CONTROL="./ams_control.sf"
WORK_DIR="./xcelium_work"

# --- Create work directory ---
mkdir -p $WORK_DIR

# --- Run xrun ---
xrun \
  -ams \
  -analogcontrol $AMS_CONTROL \
  -spectre_args "+escchars +work ahdlcmi +scsfile $NETLIST" \
  $VERILOG \
  -top SAR_1024_top \
  -access +rwc \
  -log xrun.log \
  -work $WORK_DIR
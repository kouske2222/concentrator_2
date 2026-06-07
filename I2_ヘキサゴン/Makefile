# ============================================================
# Makefile for closed-cap / open-end concentrator PO/IPO solver
# ============================================================
#
# Main usage:
#
#   make
#       Compile the Fortran solver.
#
#   make run
#       Compile and run the solver.
#
#   make visualize_bins
#       Visualize latest field_map*.bin files using visual.py.
#
#   make visualize_bins RESULT_DIR=results/260522_1530
#       Visualize field_map*.bin files in a specified result directory.
#
#   make clean
#       Remove build products only: object files, module files, executable.
#       This does NOT remove checkpoint data.
#
#   make clean_checkpoints
#       Remove IPO checkpoint directories.
#
#   make clean_results
#       Remove results directory.
#
#   make distclean
#       Remove build products, checkpoints, and results.
#
# ============================================================


# ============================================================
# Compiler settings
# ============================================================

# Fortran compiler
FC = gfortran

# Python executable.
# If .venv/bin/python exists, use it automatically.
PYTHON ?= python3
ifneq ($(wildcard .venv/bin/python),)
PYTHON := .venv/bin/python
endif

# Output directories for compilation
OUT_DIR = out
MOD_DIR = mod
BIN_DIR = bin

# Fortran compiler flags
#
# -O3:
#   Optimization.
#
# -march=native:
#   Use CPU-specific instructions on the current machine.
#
# -fopenmp:
#   Enable OpenMP parallelization.
#
# -ffree-line-length-none:
#   Allow long Fortran source lines.
#
# -J$(MOD_DIR):
#   Put .mod files into MOD_DIR.
#
# -I$(MOD_DIR):
#   Search MOD_DIR for module files.
#
FFLAGS = -O3 -march=native -fopenmp -ffree-line-length-none -J$(MOD_DIR) -I$(MOD_DIR)


# ============================================================
# Source files
# ============================================================

# Fortran source files.
# The order matters because modules must be compiled before files that use them.
SRC = mod_types.f90 \
      mod_config.f90 \
      mod_checkpoint.f90 \
      mod_geometry.f90 \
      mod_incident.f90 \
      mod_field_integrals.f90 \
      mod_currents_solver.f90 \
      mod_io.f90 \
      main_closed_cap.f90

# Object files generated from SRC
OBJS = $(SRC:%.f90=$(OUT_DIR)/%.o)

# Main executable
TARGET = $(BIN_DIR)/closed_cap_hexagonal


# ============================================================
# Checkpoint directories
# ============================================================

# Checkpoint directories used by this code.
#
# Important:
#   Do NOT use clean_checkpoints unless you really want to recompute IPO currents.
#   Checkpoints can be expensive to regenerate because humanity invented
#   repeated numerical suffering and called it simulation.
#
CHECKPOINT_DIRS = ipo_checkpoint ipo_checkpoint_closed_cap_hexagonal


# ============================================================
# Visualization settings
# ============================================================

# Python visualization script for field_map*.bin files.
#
# This script should accept:
#   --latest
#   --dir
#   --plot
#   --format
#   --dpi
#   --cmap
#   --db-floor
#   --color-percentile
#   --exclude-peak-from-colormap
#   --exclude-wall
#   --wall-radius
#   --wall-margin
#
VISUAL ?= visual.py

# Result directory control.
#
# Empty:
#   make visualize_bins
#   -> use latest directory under results/
#
# Specified:
#   make visualize_bins RESULT_DIR=results/260522_1530
#   -> use that directory
#
RESULT_DIR ?=

# Root directory searched when --latest is used.
RESULT_ROOT ?= results_hexagonal

# Plot mode for visualization.
#
# Choices expected by visual.py:
#   db      : dB plot only
#   linear  : normalized linear plot only
#   raw     : raw |E|^2 plot only
#   both    : db + linear
#   all     : db + linear + raw
#
PLOT ?= both

# Colormap.
#
# Current preference:
#   rainbow
#
CMAP ?= rainbow

# Figure output format.
#
# Choices expected by visual.py:
#   svg, png, pdf
#
FIG_FORMAT ?= png

# DPI for png/pdf.
# For svg this does not matter much, but passing it is harmless.
DPI ?= 200

# dB floor for dB plots.
#
# Example:
#   -60 means values below -60 dB are clipped to -60 dB.
#
DB_FLOOR ?= -60

# Percentile used for color scale maximum.
#
# COLOR_PERCENTILE=99:
#   The colorbar upper limit is set to the 99th percentile.
#   This helps reveal interference fringes without letting one peak dominate.
#
# COLOR_PERCENTILE=95:
#   Stronger clipping, more useful when fringes are still hidden.
#
# Empty:
#   make visualize_bins COLOR_PERCENTILE=
#   -> do not use percentile clipping.
#
COLOR_PERCENTILE ?= 99

# Exclude only the single largest value from color scale decision.
#
# 0:
#   Do not exclude the single largest point.
#
# 1:
#   Exclude the single largest point when deciding colorbar limits.
#
# Usually COLOR_PERCENTILE=99 is more useful than this,
# but this is convenient when one isolated numerical spike ruins everything.
#
EXCLUDE_PEAK ?= 0

# Exclude near-wall points when computing reference maximum and displayed peak.
#
# 0:
#   Use all points.
#
# 1:
#   Use only internal points satisfying approximately:
#     r < WALL_MARGIN * WALL_RADIUS
#
EXCLUDE_WALL ?= 0

# Wall radius used for near-wall exclusion.
# Usually this should match cfg%r_pipe.
WALL_RADIUS ?= 0.0280

# Wall margin used for near-wall exclusion.
#
# Example:
#   WALL_MARGIN=0.90 means points with r < 0.9*r_pipe are used.
#
WALL_MARGIN ?= 0.90


# ============================================================
# Phony targets
# ============================================================

.PHONY: all dirs run visualize_bins clean clean_checkpoints clean_results distclean help


# ============================================================
# Default target
# ============================================================

# make
#   Create directories and compile the main executable.
all: dirs $(TARGET)


# ============================================================
# Directory creation
# ============================================================

# make dirs
#   Create output directories required for compilation.
dirs:
	mkdir -p $(OUT_DIR) $(MOD_DIR) $(BIN_DIR)


# ============================================================
# Build rules
# ============================================================

# Link object files into the executable.
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

# Compile each Fortran source file into an object file.
#
# $<:
#   input source file
#
# $@:
#   output object file
#
$(OUT_DIR)/%.o: %.f90 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

# Explicit module dependencies.
# Without these, `make -j` can try compiling a module before its `.mod` file exists.
# Numerical suffering is bad enough without make inventing new race conditions.
$(OUT_DIR)/mod_config.o: $(OUT_DIR)/mod_types.o
$(OUT_DIR)/mod_checkpoint.o: $(OUT_DIR)/mod_types.o
$(OUT_DIR)/mod_geometry.o: $(OUT_DIR)/mod_types.o
$(OUT_DIR)/mod_incident.o: $(OUT_DIR)/mod_types.o
$(OUT_DIR)/mod_field_integrals.o: $(OUT_DIR)/mod_types.o
$(OUT_DIR)/mod_io.o: $(OUT_DIR)/mod_types.o
$(OUT_DIR)/mod_currents_solver.o: $(OUT_DIR)/mod_types.o $(OUT_DIR)/mod_checkpoint.o $(OUT_DIR)/mod_incident.o $(OUT_DIR)/mod_field_integrals.o
$(OUT_DIR)/main_closed_cap.o: $(OUT_DIR)/mod_types.o $(OUT_DIR)/mod_config.o $(OUT_DIR)/mod_geometry.o $(OUT_DIR)/mod_incident.o $(OUT_DIR)/mod_field_integrals.o $(OUT_DIR)/mod_currents_solver.o $(OUT_DIR)/mod_io.o


# ============================================================
# Run solver
# ============================================================

# make run
#   Compile the solver if necessary, then run it.
#
# TZ=Asia/Tokyo:
#   Makes timestamped output directories use Japan time.
#
run: all
	TZ=Asia/Tokyo ./$(TARGET)


# ============================================================
# Visualize binary field maps
# ============================================================

# make visualize_bins
#   Visualize latest results/*/field_map*.bin files.
#
# make visualize_bins RESULT_DIR=results/260522_1530
#   Visualize field maps in the specified directory.
#
# make visualize_bins PLOT=linear COLOR_PERCENTILE=95
#   Make linear plots with stronger color clipping to see interference fringes.
#
# make visualize_bins EXCLUDE_WALL=1
#   Compute reference value and maximum position excluding near-wall points.
#
# make visualize_bins EXCLUDE_PEAK=1
#   Exclude only the largest point from color-scale decision.
#
visualize_bins:
ifeq ($(strip $(RESULT_DIR)),)
	$(PYTHON) $(VISUAL) --latest \
		--root $(RESULT_ROOT) \
		--plot $(PLOT) \
		--format $(FIG_FORMAT) \
		--dpi $(DPI) \
		--cmap $(CMAP) \
		--db-floor $(DB_FLOOR) \
		$(if $(COLOR_PERCENTILE),--color-percentile $(COLOR_PERCENTILE),) \
		$(if $(filter 1,$(EXCLUDE_PEAK)),--exclude-peak-from-colormap,) \
		$(if $(filter 1,$(EXCLUDE_WALL)),--exclude-wall --wall-radius $(WALL_RADIUS) --wall-margin $(WALL_MARGIN),)
else
	$(PYTHON) $(VISUAL) --dir $(RESULT_DIR) \
		--plot $(PLOT) \
		--format $(FIG_FORMAT) \
		--dpi $(DPI) \
		--cmap $(CMAP) \
		--db-floor $(DB_FLOOR) \
		$(if $(COLOR_PERCENTILE),--color-percentile $(COLOR_PERCENTILE),) \
		$(if $(filter 1,$(EXCLUDE_PEAK)),--exclude-peak-from-colormap,) \
		$(if $(filter 1,$(EXCLUDE_WALL)),--exclude-wall --wall-radius $(WALL_RADIUS) --wall-margin $(WALL_MARGIN),)
endif


# ============================================================
# Cleaning targets
# ============================================================

# make clean
#   Remove only compiled files.
#
# This keeps:
#   results/
#   ipo_checkpoint/
#   ipo_checkpoint_closed_cap/
#
# Use this when you only changed source code and want to rebuild.
#
clean:
	rm -rf $(OUT_DIR) $(MOD_DIR) $(BIN_DIR)


# make clean_checkpoints
#   Remove IPO checkpoint data.
#
# Use this after changing geometry, mesh, reflection order logic,
# or anything that makes old surface-current checkpoints invalid.
#
# Do not use this casually unless you enjoy recomputing expensive things.
#
clean_checkpoints:
	rm -rf $(CHECKPOINT_DIRS)


# make clean_results
#   Remove all result directories.
#
# Use this only when old figures and field_map*.bin files are no longer needed.
#
clean_results:
	rm -rf results


# make distclean
#   Remove compiled files, checkpoints, and results.
#
# This is the nuclear option.
#
distclean: clean clean_checkpoints clean_results


# ============================================================
# Help
# ============================================================

# make help
#   Show common commands.
help:
	@echo ""
	@echo "Available targets:"
	@echo ""
	@echo "  make"
	@echo "      Compile the Fortran solver."
	@echo ""
	@echo "  make run"
	@echo "      Compile and run the solver."
	@echo ""
	@echo "  make visualize_bins"
	@echo "      Visualize latest results/*/field_map*.bin files."
	@echo ""
	@echo "  make visualize_bins RESULT_DIR=results/260522_1530"
	@echo "      Visualize field maps in a specified result directory."
	@echo ""
	@echo "  make visualize_bins PLOT=linear COLOR_PERCENTILE=95"
	@echo "      Plot linear maps with percentile color clipping."
	@echo ""
	@echo "  make visualize_bins EXCLUDE_WALL=1"
	@echo "      Exclude near-wall points when deciding reference max and peak."
	@echo ""
	@echo "  make clean"
	@echo "      Remove compiled files only."
	@echo ""
	@echo "  make clean_checkpoints"
	@echo "      Remove IPO checkpoint directories."
	@echo ""
	@echo "  make clean_results"
	@echo "      Remove results directory."
	@echo ""
	@echo "  make distclean"
	@echo "      Remove build products, checkpoints, and results."
	@echo ""
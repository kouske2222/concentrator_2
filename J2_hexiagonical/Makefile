FC := gfortran
PYTHON ?= ../.venv/bin/python
OMP_THREADS ?= $(shell nproc)
OMP_ENV := OMP_NUM_THREADS=$(OMP_THREADS) OMP_PROC_BIND=close OMP_PLACES=cores
FFLAGS ?= -O3 -march=native -fopenmp -ffree-line-length-none -Wall -Wextra -Wimplicit-interface
BUILD := build
MOD := $(BUILD)/mod
COMMON := common/mod_types.f90 common/mod_config.f90 common/mod_geometry.f90 \
          common/mod_incident.f90 common/mod_operator.f90
SOLVERS := full3d/mod_full3d.f90 c6_modal/mod_c6_modal.f90
COMPARE := $(BUILD)/compare_full_vs_c6
SWEEP := $(BUILD)/run_94ghz_sweep

.NOTPARALLEL:

.PHONY: all run run-modal run-94 preflight-94 plot plot-modes plot-all check-python install-python-deps clean debug
all: $(COMPARE) $(SWEEP)

$(COMPARE): $(COMMON) $(SOLVERS) compare_full_vs_c6.f90
	mkdir -p $(BUILD) $(MOD)
	$(FC) $(FFLAGS) -J$(MOD) -I$(MOD) $(COMMON) $(SOLVERS) compare_full_vs_c6.f90 -o $@

$(SWEEP): common/mod_types.f90 common/mod_config.f90 run_94ghz_sweep.f90
	mkdir -p $(BUILD) $(MOD)
	$(FC) $(FFLAGS) -J$(MOD) -I$(MOD) common/mod_types.f90 common/mod_config.f90 run_94ghz_sweep.f90 -o $@

run: $(COMPARE)
	$(OMP_ENV) ./$(COMPARE) config_validation_full.nml
	$(OMP_ENV) ./$(COMPARE) config_validation_modal.nml
	$(OMP_ENV) ./$(COMPARE) config_validation_post.nml

run-modal: $(COMPARE)
	$(OMP_ENV) ./$(COMPARE) config_validation_modal.nml

run-94: $(COMPARE)
	$(OMP_ENV) ./$(COMPARE) config_94ghz_full.nml
	$(OMP_ENV) ./$(COMPARE) config_94ghz_modal.nml
	$(OMP_ENV) ./$(COMPARE) config_94ghz_post.nml

preflight-94: $(SWEEP)
	./$(SWEEP) config_94ghz.nml

check-python:
	@test -x "$(PYTHON)" || { echo "ERROR: Python virtual environment not found: $(PYTHON)"; exit 1; }

install-python-deps: check-python
	$(PYTHON) -m pip install -r requirements.txt

plot: check-python
	MPLCONFIGDIR=$${MPLCONFIGDIR:-/tmp/matplotlib} $(PYTHON) plot_compare.py --result results/validation_split

plot-modes: check-python
	MPLCONFIGDIR=/tmp/matplotlib $(PYTHON) plot_modes.py --result results/validation_split

plot-all: plot plot-modes

debug: FFLAGS := -O0 -g -fopenmp -ffree-line-length-none -Wall -Wextra -Wimplicit-interface \
                  -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
debug: clean all

clean:
	rm -rf $(BUILD)

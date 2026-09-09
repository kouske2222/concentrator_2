FC := gfortran
PYTHON ?= ../.venv/bin/python
FFLAGS ?= -O3 -fopenmp -ffree-line-length-none -Wall -Wextra
COMMON := common/mod_types.f90 common/mod_config.f90 common/mod_material.f90 common/mod_geometry.f90 common/mod_incident.f90 common/mod_operator.f90
SRC := $(COMMON) c6_modal/mod_c6_modal.f90 worker.f90
.PHONY: all test preflight
all: build/worker build/libhm.so
build/worker: $(SRC)
	mkdir -p build/mod
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod $(SRC) -o $@
build/libhm.so: $(COMMON) c6_modal/mod_c6_modal.f90 acceleration/mod_bridge.f90
	mkdir -p build/hm-mod
	$(FC) $(FFLAGS) -fPIC -shared -Jbuild/hm-mod -Ibuild/hm-mod $^ -o $@
preflight:
	$(PYTHON) run.py --config config_170ghz.nml --preflight
test: all
	OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 $(PYTHON) -m unittest discover -s tests -v

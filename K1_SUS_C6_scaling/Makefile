FC := gfortran
PYTHON ?= ../.venv/bin/python
BASE := ../K1_SUS_C6_hmatrix
FFLAGS := -O3 -fopenmp -ffree-line-length-none -Wall -Wextra
SRC := $(BASE)/common/mod_types.f90 $(BASE)/common/mod_config.f90 $(BASE)/common/mod_material.f90 $(BASE)/common/mod_geometry.f90 $(BASE)/common/mod_incident.f90 $(BASE)/common/mod_operator.f90 $(BASE)/c6_modal/mod_c6_modal.f90 $(BASE)/acceleration/mod_bridge.f90 timing_bridge.f90
.PHONY: all test
all: build/libhm.so
build/libhm.so: $(SRC)
	mkdir -p build/mod
	$(FC) $(FFLAGS) -shared -fPIC -Jbuild/mod -Ibuild/mod $(SRC) -o $@
test: all
	OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 $(PYTHON) -m unittest discover -s tests -v

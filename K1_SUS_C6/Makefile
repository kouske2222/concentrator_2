FC := gfortran
PYTHON ?= python3
FFLAGS ?= -O3 -fopenmp -ffree-line-length-none -Wall -Wextra
SRC := common/mod_types.f90 common/mod_config.f90 common/mod_material.f90 common/mod_geometry.f90 common/mod_incident.f90 common/mod_operator.f90 c6_modal/mod_c6_modal.f90 worker.f90
.PHONY: all test test-python run preflight debug
all: build/worker
build/worker: $(SRC)
	mkdir -p build/mod
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod $(SRC) -o $@
preflight:
	$(PYTHON) run.py --config config_170ghz.nml --preflight
run: all
	$(PYTHON) run.py --config config_170ghz.nml --output results/170ghz
test-python:
	$(PYTHON) -m unittest discover -s tests -v
test: all test-python
	$(FC) $(FFLAGS) -Jbuild/mod -Ibuild/mod common/mod_types.f90 common/mod_config.f90 common/mod_material.f90 tests/test_surface.f90 -o build/test_surface
	./build/test_surface
	$(PYTHON) validate.py
debug:
	$(MAKE) -B FFLAGS='-O0 -g -fopenmp -ffree-line-length-none -Wall -Wextra -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace' all

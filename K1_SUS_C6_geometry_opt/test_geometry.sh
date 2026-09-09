#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")"
mkdir -p build/geometry-mod
flags="-O2 -g -fopenmp -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow"
gfortran $flags -Jbuild/geometry-mod -Ibuild/geometry-mod common/mod_types.f90 common/mod_config.f90 common/mod_material.f90 common/mod_geometry.f90 common/mod_incident.f90 common/mod_operator.f90 c6_modal/mod_c6_modal.f90 reference/mod_operator_old.f90 reference/mod_operator_exact.f90 reference/mod_modal_exact.f90 tests/test_geometry.f90 -o build/test_geometry
OMP_NUM_THREADS=2 ./build/test_geometry

make all test-python
gfortran -O3 -fopenmp -ffree-line-length-none -Jbuild/mod -Ibuild/mod common/mod_types.f90 common/mod_config.f90 common/mod_material.f90 tests/test_surface.f90 -o build/test_surface
./build/test_surface
gfortran -O3 -fopenmp -ffree-line-length-none -Jbuild/geometry-mod -Ibuild/geometry-mod common/mod_types.f90 common/mod_config.f90 common/mod_material.f90 common/mod_geometry.f90 common/mod_incident.f90 reference/mod_operator_old.f90 reference/mod_modal_old.f90 reference/worker_old.f90 -o build/worker_old
OMP_NUM_THREADS=2 ../.venv/bin/python check_geometry_run.py

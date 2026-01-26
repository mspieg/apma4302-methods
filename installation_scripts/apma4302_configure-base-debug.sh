#!/bin/bash

set -e  # Exit on any error
set -u  # Exit on undefined variable

# Script description
echo "Starting APMA 4302 PETSc configuration..."

export PETSC_DIR=$(pwd)
export PETSC_ARCH=apma4302-base-debug

./configure \
  --with-cc=mpicc \
  --with-cxx=mpicxx \
  --with-fc=mpif90 \
  --with-debugging=1 \
  --with-shared-libraries=1 \
  --COPTFLAGS="-g -O0 -march=native" \
  --CXXOPTFLAGS="-g -O0 -march=native" \
  --FOPTFLAGS="-g -O0 -march=native"


echo "Configuration completed successfully!"
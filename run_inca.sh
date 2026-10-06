#!/bin/bash

# 1. Force Linear Algebra libraries to be SERIAL (Single-threaded)
# This prevents the "hanging" issue by stopping thread explosion.
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# 2. Set Default OpenMP threads (Optional)
# If the user hasn't set OMP_NUM_THREADS, default to all available cores.
if [ -z "$OMP_NUM_THREADS" ]; then
    export OMP_NUM_THREADS=$(nproc)
fi

echo "Running INCA with OMP threads: $OMP_NUM_THREADS (BLAS forced to 1)"

# 3. Launch the actual executable
# "$@" passes all arguments (like input.inp) to the exe
./inca.exe "$@"
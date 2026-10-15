#!/bin/bash
CORES=$1
BIN=mg_cti
NAME=$2

# First load model (run0.dat) in single CPU mode to do initial refinement
# Later load time zero in parallel to better distribute the cells on each cpu
if [ ! -f ./data/$NAME.0000.h.dsb ]; then
  $BIN <<- EOF
	\$read run0
	q
	EOF
fi

# Now do the parallel computation, running run.dat
mpiexec -n $CORES $BIN <<- EOF
	\$read run
	q
	EOF



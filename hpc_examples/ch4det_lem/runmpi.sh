#!/bin/bash
CORES=$1
BIN=mg_lem


mpiexec -n $CORES $BIN <<- EOF
	\$read run
	q
	EOF





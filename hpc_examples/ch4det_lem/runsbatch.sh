#!/bin/bash
#SBATCH --ntasks=50
#SBATCH --mem-per-cpu=4096M
#SBATCH --time=7-00:00
#SBATCH --output=mgout
BIN=mg_lem

mpirun $BIN <<- EOF
        \$read run
        q
        EOF


#!/bin/bash
#SBATCH --account=def-bmaxwell
#SBATCH --mail-user=bmaxwell@uottawa.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --job-name=runsbatch.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=07-00:00
#SBATCH --output=mgout

BIN=mg_cti
NAME=$1

if [ ! -f ./data/$NAME.0000.h.dsb ]; then
  $BIN <<- EOF
	\$read run0
	q
	EOF
fi

mpirun $BIN <<- EOF
	\$read run
	q
	EOF

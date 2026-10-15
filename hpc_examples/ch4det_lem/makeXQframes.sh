#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32000M
#SBATCH --time=0-01:00
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math	# Allows simple arithmatic in bash/dash
fi


DT=$2
D_CJ=$3
NAME=$4

# get name of dat file
TMP=`grep 'DATFILE=' runsim.sh | awk -F '[=  ]' -v pat="$pattern" '{print $2}' | awk -F ' ' -v pat="$pattern" '{print $1}'`
DATFILE=$TMP.dat


# Now get boundaries
BOT=`grep 'Y ' $DATFILE | awk -F ' ' -v pat="$pattern" '{print $3}'`
LOWRBD=$(float_eval "$BOT + 0.00001")
TOP=`grep 'Y ' $DATFILE | awk -F ' ' -v pat="$pattern" '{print $4}'`
UPPERBD=$(float_eval "$TOP - 0.00001")



mg_lem <<- EOF
	paper 10.0 4.0
	norwpart read $NAME.$1
	! read pressure data along a line at the top of the domain
	mask y $UPPERBD
	plot pg export ./xq/linedata$1 p
	! Now dump 2D sootfoil data to raw xq file, adjust resolution with levs
	! be mindful of the memory requirements
	unmask *
	levs 0 4
	plot sf export ./xq/sf$1 p
	q
	EOF





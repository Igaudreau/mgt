#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0-01:00
#SBATCH --account=def-bmaxwell
#SBATCH --output=postout
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math	# Allows simple arithmatic in bash/dash
fi

NAME=$2
# get name of dat file
TMP=`grep 'DATFILE=' runsim.sh | awk -F '[=  ]' -v pat="$pattern" '{print $2}' | awk -F ' ' -v pat="$pattern" '{print $1}'`
DATFILE=$TMP.dat

# Now get boundaries
BOT=`grep 'Y ' $DATFILE | awk -F ' ' -v pat="$pattern" '{print $4}'`
LOWRBD=$(float_eval "$BOT + 0.00001")
#echo $LOWRBD
TOP=`grep 'Y ' $DATFILE | awk -F ' ' -v pat="$pattern" '{print $5}'`
UPPERBD=$(float_eval "$TOP - 0.00001")
#echo $UPPERBD


mg_cti <<- EOF
	read $NAME.$1
	! probe pressure data along bottom wall
	mask y $LOWRBD
	plot pg export ./xq/pg_bot.$1 p
	! probe pressure data along top wall
	unmask *
	mask y $UPPERBD
	plot pg export ./xq/pg_top.$1 p
	! Now dump 2D sootfoil data to raw xq file, note not currently the highest resolution, adjust with levs
	unmask *
	levs 0 4
	plot sf export ./xq/sf$1 p
	q
	EOF



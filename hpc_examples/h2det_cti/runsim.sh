#! /bin/bash
# main bash script for running mg_4s calculations
# edit file name in '_name' file
# change these parameters as necessary
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math	# Allows simple arithmatic in bash/dash
fi

NODES=1             # running on local machine with MPI directly is one node
NCPU=48			# number of cpus to use per node (distributed memory mpi)
mg_timestep=0.000002 	# time between loops
# set timestep to 0.000002 for oran case, 0.00002 for heat loss case
tstart=0.0		# time to advance before looping
LOOP_START=0		# first file to advance/read from (leave as zero always)
LOOP_END=499		# last file to advance/read from
#DATFILE=2DChannel_HeatLoss	# model dat file without the .dat extension
DATFILE=2DChannel_hotspot
NAME=det		# name of save files

# first check for existing data files
# loop from zero and set LOOP_START
# This will start the job where it left off previously
for i in $(seq -f %04g 0 $LOOP_END) ; do
  if [ -f ./data/$NAME.$i.h.dsb ]; then
    LOOP_START=$i		# first file to advance/read from
  fi
done
echo "start from loop $LOOP_START"


rm run.dat
touch run.dat
rm run0.dat
touch run0.dat
rm core*
ulimit -c unlimited
rm gmon.out
rm analysis.txt
rm jobid
cp -f $mg/mg_cti/struct.h struct.h

if [ -f mgout ]; then
  cat mgout >> mgout.log
  rm mgout
fi



# if starting from LOOP 0, load model parameters into mg and write to file
if [ "$LOOP_START" -lt "0001" ]; then
  if [ ! -f ./data/$NAME.0000.h.dsb ]; then
    mkdir data
    #load data file and write initial file, then advance with Euler method to tstart
    echo "\$read $DATFILE" >> run0.dat
    echo "write $NAME.0000" >> run0.dat
  fi
fi

# now create the run script for MG
j=0
k=0

time=0.0
LOOP_START=$(echo $LOOP_START | sed 's/^0*//')
LOOP_START=`printf %04d ${LOOP_START}`
echo "read $NAME.$LOOP_START" >> run.dat
for i in $(seq -f %04g $LOOP_START $LOOP_END) ; do
  j=$(float_eval "$i + 1")
  k=$(float_eval "$k + 1")
  time=$(float_eval "$j * $mg_timestep")
  time=$(float_eval "$tstart + $time")
  j=`printf %04d $j`
  if [ -f ./data/$NAME.$j.h.dsb ]; then
    rm -r ./data/$NAME.$j.*
  fi
  echo "target $time march 10000" >> run.dat
  echo "write $NAME.$j" >> run.dat
  # print stats every 10 steps in time
  if [ "$k" -gt "9" ]; then
    echo "stats" >> run.dat
    k=0
  fi
  # you can reload the file if you wish to defrag memory at every time interval...
  #echo "read $NAME.$j" >> run.dat
done

# check if sbatch command exists (used on cluster), otherwise run on local machine
if ! type sbatch > /dev/null; then
  echo "sbatch not found! running mpiexec directly..."
  nohup sh runmpi.sh $NCPU $NAME > mgout &
else
  # edit number of CPUS in pbs script
  sed -i '/SBATCH --nodes*/c\#SBATCH --nodes='$NODES'' runsbatch.sh
  sed -i '/SBATCH --ntasks-per-node*/c\#SBATCH --ntasks-per-node='$NCPU'' runsbatch.sh
  nohup sbatch runsbatch.sh $NAME > jobid &
fi



# Alternative syntax for MG:
#smprun -npe=$NCPU mg_4s << EOF
#nohup mpiexec -n $NCPU mg_4s < $DATFILE > mgout &
#nohup mg_4s < $DATFILE > mgout &

# OR uncomment one of these for debugging parallel run
#mpiexec -n 10 xterm -e gdb -ex run mg_4s
#gdb mg_4s


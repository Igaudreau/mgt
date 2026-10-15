#! /bin/bash
# main bash script for running mg_lem calculations
# change these parameters as necessary
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math	# Allows simple arithmatic in bash/dash
fi

NCPU=50			# number of cpus to use (distributed memory mpi)
mg_timestep=0.425	# time between loops (11us)
tstart=0.0		# time to advance before looping
LOOP_START=0		# first file to advance/read from (leave as zero always)
LOOP_END=319		# last file to advance/read from
DATFILE=detonation2D	# model dat file without the .dat extension
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
rm core*
ulimit -c unlimited
rm gmon.out
rm analysis.txt
cp -f $mg/mg_lem/struct.h struct.h	# copy this file to have a record

if [ -f mgout ]; then
  cat mgout >> mgout.log
  rm mgout
fi



# if starting from LOOP 0, load model parameters into mg and write to file
if [ "$LOOP_START" -lt "0001" ]; then
  if [ ! -f ./data/$NAME.0000.h.dsb ]; then
    # clean up from previous simulation
    #sh ./clear.sh	# too dangerous!
    mkdir data
    #load data file and write initial file, then advance with Euler method to tstart
    echo "\$read $DATFILE" >> run.dat
    echo "write $NAME.0000" >> run.dat
  fi
fi

# now create the run script for MG
j=0
k=0
time=0.0
LOOP_START=$(echo $LOOP_START | sed 's/^0*//')
LOOP_START=`printf %04d ${LOOP_START}`
echo "read $NAME.$LOOP_START" >> run.dat
echo "randomseed" >> run.dat			# needed for lem
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
  #echo "read $NAME.$j" >> run.dat
  # print stats every 10 steps in time
  if [ "$k" -gt "9" ]; then
    echo "stats" >> run.dat
    k=0
  fi

done
#echo "\n" >> run.dat


# check if srun command exists
if ! type srun > /dev/null; then
  echo "srun not found! running mpiexec directly..."
  sed -i '/SBATCH --ntasks*/c\#SBATCH --ntasks='$NCPU'' runmpi.sh
  #nohup sbatch runmpi.sh > jobid &
  nohup sh runmpi.sh $NCPU > mgout &
else
  # edit number of CPUS in pbs script
  sed -i '/SBATCH --ntasks*/c\#SBATCH --ntasks='$NCPU'' runsbatch.sh
  nohup sbatch runsbatch.sh > jobid &
fi



# Alternative syntax for MG:
#smprun -npe=$NCPU mg_lem << EOF
#nohup mpiexec -n $NCPU mg_lem < $DATFILE > mgout &
#nohup mg_lem < $DATFILE > mgout &

# OR uncomment one of these for debugging parallel run
#mpiexec -n 10 xterm -e gdb -ex run mg_lem
#gdb mg_lem


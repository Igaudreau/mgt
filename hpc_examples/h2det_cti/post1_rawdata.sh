#! /bin/bash
#####################################################################
# This is the first post-processing script
# It executes makeXQframes.sh, which extracts pressure profiles along the top
# and bottom walls, which is needed to calculate the wave speeds
# The script also extracts raw sootfoil data to xq files for each frame
# After running this scirpt, run getspeed.sh to compile the wave speed data
# Then run sootfoil.sh to compile a single high-resolution sootfoil xq file from all sootfoil data
####################################################################
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math		# Allows simple arithmatic in bash/dash
fi
LOOP_START=0		# first file to read from
#LOOP_END=200		# last file to read from
TMP=`grep 'LOOP_END=' runsim.sh | awk -F '[=  ]' -v pat="$pattern" '{print $2}' | awk -F ' ' -v pat="$pattern" '{print $1}'`
LOOP_END=$(float_eval "$TMP + 1")

NPROCS=5		# number of processes at one time
export NAME=det

# post-process simulation frames in a loop
rm slurm*.out
mkdir xq

njob=0
for i in $(seq -f %04g $LOOP_START $LOOP_END) ; do
  while float_cond $njob' > '$NPROC; do 
    sleep 1 
  done
  if [ ! -f ./xq/sf$i.xq ]; then
   if ! type sbatch > /dev/null; then
    nohup sh ./makeXQframes.sh $i $NAME &
   else
    echo $i
    sbatch ./makeXQframes.sh $i $NAME
   fi
  fi

  njob=$(float_eval "$njob + 1")
  j=$(float_eval "$NPROCS - 1")
  if float_cond $njob' > '$j ; then
    wait
    njob=0
  fi
done
wait






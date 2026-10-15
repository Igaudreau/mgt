#! /bin/bash
#####################################################################
# This is the second post-processing script
# It executes makePNGframes.sh, which plots 2D figures of
# temperature, density, pressure, density gradient, schlieren
# After running this script, create a movie by executing makemovie.sh
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
D_CJ=`grep 'Dcj =' mgout | awk -F ' ' -v pat='$pattern' '{print $3}'`
mg_timestep=`grep 'mg_timestep=' runsim.sh | awk -F '[=  ]' -v pat="$pattern" '{print $2}' | awk -F ' ' -v pat="$pattern" '{print $1}'`
echo mg_timestep $mg_timestep
echo D_CJ $D_CJ


# post-process simulation frames in a loop
rm slurm*.out
mkdir png

njob=0
for i in $(seq -f %04g $LOOP_START $LOOP_END) ; do
  while float_cond $njob' > '$NPROC; do 
    sleep 1 
  done
  if [ ! -f ./png/sch$i.png ]; then
   if ! type sbatch > /dev/null; then
    nohup sh ./makePNGframes.sh $i $NAME $D_CJ $mg_timestep &		# call mg_post.sh
   else
    echo $i
    sbatch ./makePNGframes.sh $i $NAME $D_CJ $mg_timestep
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



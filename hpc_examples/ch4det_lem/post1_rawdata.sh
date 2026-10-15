#! /bin/bash
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math		# Allows simple arithmatic in bash/dash
fi

NPROCS=8		# number of processes at one time
D_CJ=6.3		# analytical solution to this problem

# get variables from runsim.sh
LOOP_START=0		# first file to read from
TMP=`grep 'LOOP_END=' runsim.sh | awk -F '[=  ]' -v pat="$pattern" '{print $2}' | awk -F ' ' -v pat="$pattern" '{print $1}'`
LOOP_END=$(float_eval "$TMP + 1")
echo LOOP_END $LOOP_END
mg_timestep=`grep 'mg_timestep=' runsim.sh | awk -F '[=  ]' -v pat="$pattern" '{print $2}' | awk -F ' ' -v pat="$pattern" '{print $1}'`
echo mg_timestep $mg_timestep
export NAME=`grep 'NAME=' runsim.sh | awk -F '[=  ]' -v pat="$pattern" '{print $2}' | awk -F ' ' -v pat="$pattern" '{print $1}'`
echo NAME $NAME

# post-process simulation frames in a loop
rm slurm*.out
mkdir xq

njob=0
for i in $(seq -f %04g $LOOP_START $LOOP_END) ; do
  while float_cond $njob' > '$NPROC; do 
    sleep 1 
  done
  if [ ! -f ./xq/u3$i.xq ]; then
   if ! type sbatch > /dev/null; then
    nohup sh ./makeXQframes.sh $i $mg_timestep $D_CJ $NAME &		# call mg_post.sh
   else
    sbatch ./makeXQframes.sh $i $mg_timestep $D_CJ $NAME
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






#Useful on HPC to check status of jobs


#! /bin/bash

HOME=/home/smiri

cd $HOME/scratch
for i in ./5cylinders/*/; do
  echo ---------
  echo $i:
  #cd $i
  cat ${i}jobid
  #cd data
  for j in ${i}data/det.*.h.dsb ; do
    k=$j
  done
  echo $k
  for j in ${i}core* ; do
    k=$j
  done
  echo $k
  #cd ..
  #cd ..
done

#squeue -u bmaxwell -o "%.18i %C  %D %.8j %.8u %.2t %.10M %r"
sq
quota -u
~             
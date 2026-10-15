#!/bin/bash

# this script will copy the clean mg_g to a new folder and setup the environment variables etc..
#

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "Syntax is 'sh ./removemg.sh name' where name is the version name of mg to delete (i.e. mg_name)"
elif [ -d "./mg_$1" ]; then
    # do not remove mg_g
    if [ $1 != "g" ]; then
      read -p "Are you sure you want to remove mg_$1 completely? (y/n)" REPLY
      #read -p "Are you sure you want to remove mg_$1 completely? (y/n)" -n 1 -r 
      echo    # move to a new line
      #if [[ $REPLY =~ ^[Yy]$ ]]; then
      if [ "$REPLY" = "y" ]; then
         # do dangerous stuff
         echo "removing mg_$1..."
         # delete files and directories
         rm -r ./mg_$1/*
         rmdir mg_$1
         # delete environment variables
         sed -i 's#$mg/mg_'$1'/bin ##g' _cshrc
         sed -i '\#mg_'$1'#d' _cshrc
         sed -i 's#$mg/mg_'$1'/bin:##g' _bashrc
         sed -i '\#mg_'$1'#d' _bashrc
         #sed -i '\#mg_'$1'#d' installbash.sh
         sed -i '\#if (-d mg_'$1') then#,+7 d' Clearall
         sed -i '\#if (-f bin/mg_'$1') then#,+2 d' Clearall
         sed -i '\#cd $mg/mg_'$1'#,+1 d' Clearall
         echo "mg_$1 deleted.  Please remember to remove this version from the git repo."

      fi
    else
      echo "error: cannot remove mg_g (by default)."
    fi
else
  echo "error: mg_$1 not found!"

fi







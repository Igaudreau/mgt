#! /bin/bash

# this script will copy the mg_g program to a new folder and setup the environment variables etc..
#

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "Syntax is 'sh ./newmg.sh name' where name is the version name of mg to create (i.e. mg_name)"
elif [ -d "./mg_$1" ]; then
    echo "mg_$1 already exists!"
else
  echo "creating mg_$1..."
  # first copy mg_g
  mkdir mg_$1
  cp -r mg_g/* mg_$1/

  cp -f  _bashrc  _bashrc.bak
  cp -f  _cshrc  _cshrc.bak
  # add new version to environment variables
#  sed -i '\#mkdir $mg/mg_g/data#a mkdir $mg/mg_'$1'/data' install.sh
#  sed -i '\#mkdir $mg/mg_g/data#a mkdir $mg/mg_'$1'/bin' install.sh
  sed -i '\#  setenv mg_g_data $mg/mg_g/data#a\  setenv mg_'$1'_data $mg/mg_'$1'/data' _cshrc
  sed -i 's#set path=($path#set path=($path $mg/mg_'$1'/bin#g' _cshrc

  #sed -i '\#export mg_g_data=$mg/mg_g/data#aexport mg_'$1'_data=$mg/mg_'$1'/data' _bashrc
  sed -i '1i export mg_'$1'_data=./data' _bashrc
  sed -i 's#PATH=#PATH=$mg/mg_'$1'/bin:#g' _bashrc
  #source _bashrc

  sed -i '\#cd mg_g#a\  cd mg_'$1'' Clearall
  sed -i '\#cd mg_g#a if (-d mg_'$1') then' Clearall
  sed -i '\#cd mg_g#a endif' Clearall
  sed -i '\#cd mg_g#a\  cd ..' Clearall
  sed -i '\#cd mg_g#a\  endif' Clearall
  sed -i '\#cd mg_g#a\    ./Clearsmp' Clearall
  sed -i '\#cd mg_g#a\  if (-d .smp) then' Clearall
  sed -i '\#cd mg_g#a\  ./Clear' Clearall

  sed -i '\#rm bin/mg_g#a\  rm bin/mg_'$1'' Clearall
  sed -i '\#rm bin/mg_g#a if (-f bin/mg_'$1') then' Clearall
  sed -i '\#rm bin/mg_g#a endif' Clearall
  sed -i '\#cd $mg/mg_g#a cd $mg/mg_'$1'' Clearall
  sed -i '\#cd $mg/mg_g#a rm -r *~' Clearall
  mkdir $mg/mg_$1/data
  mkdir $mg/mg_$1/bin
  sed -i 's#mg_g#mg_'$1'#g' $mg/mg_$1/makemg
  sed -i 's#mg_g#mg_'$1'#g' $mg/mg_$1/Clear
  sed -i 's#mg_g#mg_'$1'#g' $mg/mg_$1/usr.c
  # create uppercase version
  UCAS1=$( echo "$1" | tr -s  '[:lower:]'  '[:upper:]' )
  sed -i 's#MG_G#MG_'$UCAS1'#g' $mg/mg_$1/usr.c
  #source ~/.bashrc

  echo "mg_$1 created.  Please remember to push the new version to git repo. You might have to re-source .bashrc in your home folder"

fi







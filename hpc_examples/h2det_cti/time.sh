#! /bin/bash
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math		# Allows simple arithmatic in bash/dash
fi

# reads mgout to figure out how much total time has elapsed
#rm timeout
#echo "Simulation time (s) = "
grep 'that took' mgout | awk -v pat="$pattern" '{print $3}' > timeout
TIME_SEC=`awk '{ sum+=$1} END {print sum}' < timeout`
TIME_MIN=$(float_eval "$TIME_SEC / 60")
TIME_HR=$(float_eval "$TIME_MIN / 60")
echo "Total simulation time = " $TIME_SEC " s /  " $TIME_MIN " min /  " $TIME_HR " hrs"



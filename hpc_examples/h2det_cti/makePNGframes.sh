#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00
#SBATCH --account=def-bmaxwell
#SBATCH --output=postout
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math	# Allows simple arithmatic in bash/dash
fi


DPI=600			# resolution


Frame=$1
NAME=$2
D_CJ=1587.01
#D_CJ=$3
DT=$4
TIME=$(float_eval "$DT * $Frame")
LOC=$(float_eval "$D_CJ * $TIME")
XL=$(float_eval "$LOC - 0.12") #Left boundry
XR=$(float_eval "$LOC + 0.06") #Right boundry


mg_cti <<- EOF
	paper 10.0 4.0
	read $NAME.$1
	! zoom in as needed
	xr $XL $XR
	! Limit size of file by only plotting coarser levs
	!levs 0 5

	grey_rain
	open a$1 ps_rgb
	nocontour
	qr 0.084 0.6
	plot rho p
	close a$1

	nogrey_rain
	open b$1 ps_rgb
	nocontour
	qr 300 3300
	plot te p
	close b$1

	grey_rain
	open c$1 ps_rgb
	nocontour
	qr 0 100
	plot rgrad p
	close c$1

	grey_rain
	open d$1 ps_rgb
	nocontour
	qr -100 100
	plot sch p
	close d$1

	q
	EOF



plot a$1 rho$1
plot b$1 te$1
plot c$1 rgrad$1
plot d$1 sch$1


for var in rho te rgrad sch ; do
  # 1. Render at 300 DPI (no need to go higher if the gaps are in the vector)
  gs -q -dSAFER -dBATCH -dNOPAUSE -dNOPROMPT \
     -sDEVICE=pngalpha -r$DPI \
     -o ./png/tmp$1.png $var$1.ps
  # 2. Heal transparent lines by dilating the alpha channel (IM6 compatible)
  convert ./png/tmp$1.png \
      -channel A -morphology Dilate Square:1 +channel \
      -background white -flatten -trim +repage ./png/$var$1.png
  rm ./png/tmp$1.png
done


rm a$1.ps
rm b$1.ps
rm c$1.ps
rm d$1.ps
rm ./png/tmp$1.png
rm rho$1.ps
rm te$1.ps
rm rgrad$1.ps
rm sch$1.ps




#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32000M
#SBATCH --time=0-01:00
if [ "$BASH_VERSION" = '' ]; then
  . ./_math		# Allows simple arithmatic in bash/dash
else
  source ./_math	# Allows simple arithmatic in bash/dash
fi


DT=$2
D_CJ=$3
NAME=$4
TIME=$(float_eval "$DT * $1")
LOC=$(float_eval "$D_CJ * $TIME")
XL=$(float_eval "$LOC + 0")
XR=$(float_eval "$LOC + 50")
#AREA=`grep 'Z ' detonation3D.dat | awk -v pat="$pattern" '{print $4}'`


DPI=600			# resolution


#srun mg_lem <<- EOF
mg_lem <<- EOF
	paper 10.0 4.0
	xr $XL $XR
	norwpart read $NAME.$1
	levs 0 5

	! density
	grey_rain
	open a$1 ps_rgb
	nocontour
	qr 0 10
	plot rho p
	close a$1

	! temperature
	nogrey_rain
	open b$1 ps_rgb
	nocontour
	qr 0.85 12
	plot te p
	close b$1

	! schlieren image
	grey_rain
	open c$1 ps_rgb
	nocontour
	qr -5 5
	plot sch p
	close c$1

	! density gradient
	grey_rain
	open d$1 ps_rgb
	nocontour
	qr 0 5
	plot rgrad p
	close d$1
	q
	EOF

plot a$1 rho$1
plot b$1 te$1
plot c$1 sch$1
plot d$1 rgrad$1


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




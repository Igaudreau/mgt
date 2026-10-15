#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32000M
#SBATCH --time=0-01:00
#SBATCH --account=def-bmaxwell
#SBATCH --output=postout
if [ "$BASH_VERSION" = '' ]; then
  . ./_math             # Allows simple arithmatic in bash/dash
else
  source ./_math                # Allows simple arithmatic in bash/dash
fi

DPI=600

mkdir png
rm sootfoil_partial.png
rm sootfoil_partial.ps
rm sootfoil_full.png
rm sootfoil_full.ps
rm sootfoil.xq

TMP=`grep 'LOOP_END=' runsim.sh | awk -F '[=  ]' -v pat="$pattern" '{print $2}' | awk -F ' ' -v pat="$pattern" '{print $1}'`
LOOP_END=$(float_eval "$TMP + 1")
python3 sootfoil.py --l $LOOP_END


mg_cti <<- EOF
	paper 10.0 4.0
	open d ps_rgb
	nocontour
	qr 130000 600000
	xr -1 1.55 readxq sootfoil
	close d
	open b ps_rgb
	nocontour
	qr 130000 600000
	xr 1.35 1.55 readxq sootfoil
	close b
	q
	EOF



plot b sootfoil_partial
plot d sootfoil_full

# sootfoil
for var in sootfoil_full sootfoil_partial ; do
  # 1. Render at 300 DPI (no need to go higher if the gaps are in the vector)
  gs -q -dSAFER -dBATCH -dNOPAUSE -dNOPROMPT \
     -sDEVICE=pngalpha -r$DPI \
     -o ./png/tmp.png $var.ps
  # 2. Heal transparent lines by dilating the alpha channel (IM6 compatible)
  convert ./png/tmp.png \
      -channel A -morphology Dilate Square:1 +channel \
      -background white -flatten -trim +repage ./png/$var.png
  rm ./png/tmp.png
done

rm sootfoil_full.ps
rm sootfoil_partial.ps
rm b.ps
rm d.ps



set terminal postscript eps enhanced color font "Times-Roman,16"
set datafile separator " "
set output "shockTube_u.eps"
set title "density profile, t=0.0008s"
set size 0.5,0.5
set yrange [0:5e8]
set xrange [-5:5]
set xlabel "x (km)"
set ylabel "u (cm/s)"
set format y "{%g}"
set format x "{%g}"
set key top
#set key outside
plot "shocktube_results.dat" using ($1/100000):3 with lines ls 7 lw 3 title 'Exact',\
"sod_u.xq" every 2::3 using ($1/100000):2 with points ls 6 title 'Roe'

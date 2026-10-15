set terminal postscript eps enhanced monochrome
set datafile separator " "
set output "Dspeed.eps"
#set title "Detonation Speed"
set size 0.8,0.8
set yrange [0.5:10]
set y2range [178.18:3563.6]
set xrange [0:0.0035]
set xlabel "Simulation time [s]"
set ylabel "Detonation Mach number ({/Times-Italic M}_D)"
set y2label "Detonation velocity [m/s]"
set ytics nomirror
set y2tics nomirror
set key at 0.0032,3
set style line 1 lt 2 lw 3 lc -1 dt 1
set style line 2 lt 3 lw 0.5 pt 6 ps 0.8 lc 3
set style line 3 lt 3 lw 5 lc 1
set style line 4 lt 4 lw 3 lc 9
set style line 5 lt 1 lw 5 lc 3 dt 6
set style line 6 lt 4 lw 3 lc 9 dt 3
D(x) = 6.35	# old simulations
S(x) = 6.3
C(x) = 5.19
B(x) = 5.53
h(x) = mean_y
fit [50:130] h(x) "shockspeed.txt" using 2:1 via mean_y
k(x) = ( 500 < x ) ? h(x) : 1/0 
plot S(x) with lines ls 1 title '{/Times-Italic M}_{CJ} = 6.30',\
"shockspeed.txt" using ($2*0.000027164):1 with linespoints ls 2 title '{/Times-Italic M}_{D}',\
h(x) with lines ls 5 title sprintf("Mean: {/Times-Italic M}_D=%.2f", mean_y),\
B(x) with lines ls 3 title '{/Times-Italic M}_{D,avg,exp},Kiyanda and Higgins} = 5.53',\
D(x) with lines ls 6 title '{/Times-Italic M}_{D,avg,LES},Maxwell et al.} = 6.35',\
C(x) with lines ls 4 title '{/Times-Italic M}_{D,avg,exp},Maxwell et al.} = 5.19'
#g(x) with lines ls 3 title 'running mean'


set terminal postscript eps enhanced color font "Times-Roman,16"
set datafile separator " "
set output "gamezo_znd.eps"
set title "ZND structure for 0.5C+0.5O"
set size 0.8,0.8
set yrange [0:1]
set xrange [1E-3:*]
set y2range[0:10]
set y2tics 2
set xlabel "{/Times-Italic x} (cm)"
set ylabel "{/Times-Italic X}, {/Times-Italic M}"
set y2label "{/Times-Italic T} (10^9 K)"
set format y "{%g}"
set format x "{%g}"
set logscale x
set key outside below
#set key samplen 2
plot "./gamezo_znd_xa_c12_x.xq" every ::3 using 1:2 with lines ls 8 linecolor 6 title 'C^{12}',\
"./gamezo_znd_xa_o16_x.xq" every ::3 using 1:2 with lines ls 2 lw 3 linecolor 1 title 'O^{16}',\
"./gamezo_znd_xa_si28_x.xq" every ::3 using 1:2 with lines ls 1 lw 3 linecolor 2 title 'Si^{28}',\
"./gamezo_znd_xa_ni56_x.xq" every ::3 using 1:2 with lines ls 3 lw 3 linecolor 4 title 'Ni^{56}',\
"./gamezo_znd_mach_x.xq" every ::3 using 1:2 with lines ls 4 lw 3 linecolor -1 dt 3 title '{/Times-Italic M}',\
"./gamezo_znd_temp_x.xq" every ::3 using 1:($2/10e9) with lines ls 4 lw 3 linecolor 7 dt 4 title '{/Times-Italic T}'



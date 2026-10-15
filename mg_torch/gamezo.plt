set terminal postscript eps enhanced color font "Times-Roman,16"
set datafile separator " "
set output "gamezo.eps"
set title "Evolution for 0.5C+0.5O"
set size 0.8,0.8
set yrange [0:1]
set xrange [*:*]
set y2range[0:10]
set y2tics 2
set xlabel "{/Times-Italic t} (s)"
set ylabel "{/Times-Italic X}"
set y2label "{/Times-Italic T} (10^9 K)"
set format y "{%g}"
set format x "{%g}"
set logscale x
set key left
#set key samplen 2
plot "./reactor_xa_c12_t.xq" every ::3 using 1:2 with lines ls 8 linecolor 6 title 'C^{12}',\
"./reactor_xa_o16_t.xq" every ::3 using 1:2 with lines ls 2 lw 3 linecolor 1 title 'O^{16}',\
"./reactor_xa_si28_t.xq" every ::3 using 1:2 with lines ls 1 lw 3 linecolor 2 title 'Si^{28}',\
"./reactor_xa_ni56_t.xq" every ::3 using 1:2 with lines ls 3 lw 3 linecolor 4 title 'Ni^{56}',\
"./reactor_temp_t.xq" every ::3 using 1:($2/10e9) with lines ls 4 lw 3 linecolor 7 dt 4 title '{/Times-Italic T}'


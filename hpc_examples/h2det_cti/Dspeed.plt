set terminal postscript eps enhanced monochrome font "Times,14"
set datafile separator " "
set output "Dspeed.eps"
set title "Wave speed along bottom wall"
set size 0.6,0.6
set yrange [0:4000]
set y2range [0:11.224604333]
set xrange [0:*]
set xlabel "Wave position [m]"
set ylabel "Wave velocity [m/s]"
set y2label "Mach number"
set ytics nomirror
set y2tics nomirror
set y2tics 1
set key below
set style line 1 lt 2 lw 3 lc -1 dt 1
set style line 2 lt 3 lw 0.25 pt 6 ps 0.5 lc 3
set style line 3 lt 3 lw 5 lc 1
set style line 4 lt 4 lw 3 lc 9
set style line 5 lt 1 lw 2 lc 3 dt 1
set style line 6 lt 4 lw 3 lc 9 dt 3
set style line 9 lt 9 lw 3 lc 7 dt 1
##### CODE FOR AVERAGE ###########
h(x) = mean_y
fit [0.2:*] h(x) "shockspeed_bot.txt" using 2:1 via mean_y
k(x) = ( 0.2 < x ) ? h(x) : 1/0
##### END CODE FOR AVERAGE ###########
##### CODE FOR MOVING AVERAGE ###########
n = 50
do for [i=1:n] {
    eval(sprintf("back_n%d=0", i))
}
shift = "("
do for [i=n:2:-1] {
    shift = sprintf("%sback_n%d = back_n%d, ", shift, i, i-1)
}
shift = shift."back_n1 = x)"
sum = "(back_n1"
do for [i=2:n] {
    sum = sprintf("%s+back_n%d", sum, i)
}
sum = sum.")"
samples(x) = $0 > (n-1) ? n : ($0+1)
avg_n(x)   = (shift_n(x), @sum/samples($0))
shift_n(x) = @shift
###### END CODE FOR MOVING AVERAGE #######
# CJ-Detonation Speed (get from mgout)
D(x) = real(system("grep 'Dcj =' mgout | awk -F ' ' -v pat='$pattern' '{print $3}'"))
plot D(x) with lines ls 1 title 'Detonation {/Times-Italic D}_{CJ}',\
"shockspeed_bot.txt" using 2:1 with linespoints ls 2 title '{/Times-Italic D}',\
"shockspeed_bot.txt" using 2:(avg_n($1)) with lines ls 5 title 'Moving Avg',\
k(x) with lines ls 3 title sprintf("Mean: {/Times-Italic D}=%.2f", mean_y)


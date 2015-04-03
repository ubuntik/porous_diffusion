#!/opt/local/bin/gnuplot
#

# set terminal svg size 350,262 fname 'Verdana' fsize 10
# set output 'introduction.svg'

reset
# png
set terminal pngcairo size 1024,768 enhanced font 'Verdana,14'
set output 'sample_9.png'
# svg
#set terminal svg size 1024,768 fname 'Verdana, Helvetica, Arial, sans-serif' \
#fsize '9' rounded dashed
#set output 'plot.svg'

# define axis
# remove border on top and right and set color to gray
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror
# define grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12

# color definitions
set style line 1 lc rgb 'red' pt 1 ps 0.5 lt 1 lw 3 # --- red
set style line 2 lc rgb 'orange' pt 1 ps 0.25 lt 1 lw 2 # --- orange
set style line 3 lc rgb 'yellow' pt 1 ps 0.25 lt 1 lw 2 # --- yellow
set style line 4 lc rgb 'green' pt 1 ps 0.25 lt 1 lw 2 # --- green
set style line 5 lc rgb 'grey' pt 1 ps 0.25 lt 1 lw 2 # --- dodgerblue
set style line 6 lc rgb 'blue' pt 1 ps 0.25 lt 1 lw 2 # --- blue
set style line 7 lc rgb 'purple' pt 1 ps 0.25 lt 1 lw 2 # --- purple
set style line 8 lc rgb 'black' pt 1 ps 0.25 lt 1 lw 2 # --- indigo

set key top right

set title 'Sample #9'
set xlabel 'x'
#set ylabel 'p(t)'
set ylabel "p(x,{/Symbol z})"
set xrange [0:1000]
set yrange [0:1]


plot 'plot.dat' u 1:2 t 'p1' w lp ls 1, \
	""	u 1:3 t 'p2' w lp ls 2, \
	""	u 1:4 t 'p3' w lp ls 3, \
	""	u 1:5 t 'p4' w lp ls 4, \
	""	u 1:6 t 'p5' w lp ls 5, \
	""	u 1:7 t 'p6' w lp ls 6, \
	""	u 1:8 t 'p7' w lp ls 7, \
	""	u 1:9 t 'p8' w lp ls 8

set output "transport_ssl.eps"
set encoding iso_8859_1
# letter is 21.59 x 27.94 cm
# a4 is 21 x 29.7 cm
# general format is 1 cm for bottom and top, 0.25 at the
# top and 0.75 at the bottom, then divide spacing
# on as many figures as is needed, i.e. 2.2 cm high each
set terminal postscript eps enhanced color size 5cm,7.6cm
set nokey
set style line 10 lc rgb "#000000" lt -1 lw 1
set style line 1 lc rgb "#000000" lt -1 lw 2
set style line 2 lc rgb "#0000FF" lt -1 lw 2
set style line 3 lc rgb "#228B22" lt -1 lw 2
set style line 4 lc rgb "#8A2BE2" lt -1 lw 2
set style line 5 lc rgb "#000000" lt -1 lw 4 dt 2
set tics scale 1
set border ls 10
set xtics font ",12" offset 0, -0.5
set ytics font ",12"
set format x "%1.1f"
set lmargin at screen 0.3
set multiplot layout 3, 1
set title ""
set tmargin at screen 0.967; set bmargin at screen 0.678;
set ylabel "{/Symbol s} [S/m]" font ",12" offset -1.8, 0
set xrange [-0.4:0.4]
set logscale y
set format y "10^{%T}"
set yrange [5e-1:1e8]
set ytics 1e1,1e2,1e7
unset mytics
#set key samplen 2
#set key bottom left font ",12"
set xtics format ""
set xtics -0.4, 0.2, 0.4
plot 'sigma' index 0 u 1:5 with linespoints ls 2 title "100 K", '' index 2 u 1:5 with linespoints ls 3 title "300 K",  '' index 6 u 1:5 with linespoints ls 4 title "700 K"
set title ""
set tmargin at screen 0.678; set bmargin at screen 0.388
unset logscale y
unset format y
set ylabel "S [{/Symbol m}V/m]" font ",12" offset 0.1, 0
set xrange [-0.4:0.4]
set yrange [-750:750]
unset key
set ytics  -500, 250, 500
set xtics -0.4, 0.2, 0.4
set xtics format ""
plot 'seebeck' index 0 u 1:5 with linespoints ls 2, '' index 2 u 1:5 with linespoints ls 3,  '' index 6 u 1:5 with linespoints ls 4
set title ""
set tmargin at screen 0.388; set bmargin at screen 0.099
set xlabel "Chemical potential [eV]" font ",12"
set ylabel "L [10^{-8} V^2/K^2]" font ",12" offset -1, 0
set xrange [-0.4:0.4]
set yrange [1.2:4.1]
set ytics 1.4, 0.5, 3.9
set format y "%1.1f"
set xtics center offset 0 0 -0.4, 0.2, 0.4
set xtics format "%1.1f"
plot 'lorenz' index 0 u 1:5 with linespoints ls 2, '' index 2 u 1:5 with linespoints ls 3,  '' index 6 u 1:5 with linespoints ls 4

set output "transport_nc_bands.eps"
set encoding iso_8859_1
# letter is 21.59 x 27.94 cm
# a4 is 21 x 29.7 cm
# general format is 1 cm for bottom and top, 0.25 at the
# top and 0.75 at the bottom, then divide spacing
# on as many figures as is needed, i.e. 2.2 cm high each
set terminal postscript eps enhanced color size 5cm,5.4cm
set style line 10 lc rgb "#000000" lt -1 lw 1
set style line 1 lc rgb "#000000" lt -1 lw 2
set style line 2 lc rgb "#0000FF" lt -1 lw 2
set style line 3 lc rgb "#228B22" lt -1 lw 2
set style line 4 lc rgb "#8A2BE2" lt -1 lw 2
set style line 5 lc rgb "#0000FF" lt 0 lw 4
set style line 6 lc rgb "#228B22" lt 0 lw 4
set style line 7 lc rgb "#8A2BE2" lt 0 lw 4
set tics scale 1
set border ls 10
set xtics font ",12" offset 0, -0.5
set ytics font ",12"
set format x "%1.1f"
set lmargin at screen 0.3
set multiplot layout 2, 1
set title ""
set tmargin at screen 0.954; set bmargin at screen 0.546
set ylabel "n [10^{21} cm^{-1}]" font ",12" offset 0, 0
set logscale y
set xrange [-0.4:0.4]
set yrange [1e-17:1e2]
set ytics left offset -3.5, 0 1e-18, 1e3, 1e2
set format y "10^{%T}"
set xtics -0.4, 0.2, 0.4
set xtics format ""
plot 'sigma' index 0 u 1:3 with linespoints ls 2 notitle, '' index 2 u 1:3 with linespoints ls 3 notitle, '' index 6 u 1:3 with linespoints ls 4 notitle, '' index 0 u 1:4 with lines ls 5 notitle, '' index 2 u 1:4 with lines ls 6 notitle, '' index 6 u 1:4 with lines ls 7 notitle
set title ""
unset logscale
set format y "%1.1f"
set tmargin at screen 0.546; set bmargin at screen 0.139
set xlabel "Chemical potential [eV]" font ",12"
set ylabel "k-point [{\305}^{-1}]" font ",12" offset 0.3, 0
set xrange [-0.4:0.4]
set yrange [-0.5:0.5]
set ytics right offset 0 0 -0.3, 0.3, 0.3
set xtics center offset 0 0 -0.4, 0.2, 0.4
set xtics format "%1.1f"
set samples 100000
plot sqrt(1.0*(x-0.2)/3.80998194577014799921) with linespoints ls 1 notitle, -sqrt(1.0*(x-0.2)/3.80998194577014799921) with linespoints ls 1 notitle, sqrt(-0.1*(x-0.0)/3.80998194577014799921) with linespoints ls 1 notitle, -sqrt(-0.1*(x-0.0)/3.80998194577014799921) with linespoints ls 1 notitle, sqrt(-10*(x+0.3)/3.80998194577014799921) with linespoints ls 1 notitle, -sqrt(-10*(x+0.3)/3.80998194577014799921) with linespoints ls 1 notitle

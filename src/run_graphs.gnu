set terminal png font "arial,10" fontscale 1.0
set xlabel 'delta r Steps'
set ylabel 'Pair Correlation'
set grid
set key off
set title 'Pair Correlation'
set output 'pair_corre.png'
plot 'pair_corre_data.txt' with linespoints lc rgb "blue"

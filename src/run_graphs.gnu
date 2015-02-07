set terminal png font "arial,10" fontscale 1.0
set xlabel 'Time'
set ylabel 'Energy'
set grid
set key off

set title 'Potential Energy'
set output 'ener_pot.png'
plot 'ener_pot_data.txt' with linespoints lc rgb "blue"

set title 'Kinetic Energy'
set output 'ener_kin.png'
plot 'ener_kin_data.txt' with linespoints lc rgb "blue"

set title 'Total Energy'
set output 'ener_tot.png'
plot 'ener_tot_data.txt' with linespoints lc rgb "blue"

set title 'Final Temperature'
set output 'temp_final.png'
plot 'temp_final_data.txt' with linespoints lc rgb "blue"

set title 'Pair Correlation'
set output 'pair_corre.png'
plot 'pair_corre_data.txt' with linespoints lc rgb "blue"

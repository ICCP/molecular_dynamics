set terminal png font "arial,10" fontscale 1.0
set xlabel 'Time'
set ylabel 'Energy'
set grid
set key off

set title 'Potential Energy'
set output 'ener_pot.png'
plot 'ener_pot_data.txt' with linespoints ls 6 lc 7

set title 'Kinetic Energy'
set output 'ener_kin.png'
plot 'ener_kin_data.txt' with linespoints ls 6 lc 7

set title 'Total Energy'
set output 'ener_tot.png'
plot 'ener_tot_data.txt' with linespoints ls 6 lc 7

set title 'Temperature'
set output 'temperature.png'
plot 'temp_data.txt' with linespoints ls 6 lc 7
set terminal png font "arial,10" fontscale 1.0
set xlabel 'Time (x0.001 sec)'
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

set ylabel 'Temperature'
set title 'Temperature'
set output 'temperature.png'
plot 'temp_data.txt' with linespoints lc rgb "blue"
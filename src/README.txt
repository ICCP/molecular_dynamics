Compile Fortran code:
f95 -o mdfort.exe mdfortran.f

1) Please provide input to md.py
   Please ENSURE L = # of lattice is same in all three codes
2) Execute md.py to generate initial lattice
3) Execute mdfort.exe for molecular dynamics calculations
4) Again execute md.py to generate all the plots
5) Execute rdf.py for radial distribution function plot
rm *.o
rm *.mod
rm *.out
rm data/*
f95 -c initial_conditions.f95
f95 -c verlet_algorithm.f95
f95 -c main.f95
f95 verlet_algorithm.o initial_conditions.o main.o
./a.out


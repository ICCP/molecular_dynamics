clear all; close all; clc;
load ../energy.out;

NT = 500;

time = zeros(NT,1);
tot_energy = zeros(NT,1);
pe = zeros(NT,1);
ke = zeros(NT,1);

for ii = 1:NT
    time(ii,1) = energy(ii,1);
    tot_energy(ii,1) = energy(ii,2);
    pe(ii,1) = energy(ii,3);
    ke(ii,1) = energy(ii,4);
end 

S = std(tot_energy(250:500,1));

plot(time(:,1),ke(:,1));
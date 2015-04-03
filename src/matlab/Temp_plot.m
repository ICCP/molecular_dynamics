clear all; close all; clc;
load ../temperature.out;
%load ../pressure.out
NT = 500;

time = zeros(NT,1);
Temp = zeros(NT,1);

for ii = 1:NT
    time(ii,1) = temperature(ii,1);
    Temp(ii,1) = temperature(ii,2);
    
end

plot(Temp(:,1));
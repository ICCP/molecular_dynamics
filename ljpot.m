clc; clear all; close all;

R = linspace(0.5,2,2001);
ep = 1e-2;
sig = 1;

vLJ = 4 * ep * ((sig./R).^12 - (sig./R).^6);

figure;
plot(R,vLJ)
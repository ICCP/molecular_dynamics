% Plot the pair correlation function
clc;
clear all;
close all;

rr = load('fort.77');
rr = sort(rr);
maxr = max(rr);
dr = 0.1;
Nr = ceil(maxr/dr);
rdist = 0:dr:Nr*dr;
boxes = zeros(Nr,1);
for i = 1:Nr
    tmp1 = rr > rdist(i);
    tmp2 = rr <= rdist(i+1);
    cnt = sum(tmp1.*tmp2);
    boxes(i) = cnt;
end

figure;
plot(rdist(2:end),boxes);
grid on;
xlabel('Distance');
ylabel('Frequency');
title('Pair correlation function, T = 15');
clear all; close all; clc;
load ../raddist.out;
NT = 500;
Bins = 100
Dist = zeros(NT,Bins);
plotbox = zeros(1,Bins);
for ii = 1:NT
    for jj = 1:Bins
        Dist(ii,jj) = raddist((ii-1)*Bins+jj,2);
        if ii > 250
            plotbox(1,jj) = plotbox(1,jj) + raddist((ii-1)*Bins+jj,2)/250;
        end
    end
end

neighbor_dist = 2^(1/6);
boxlength = 2/2^.5;
bounds = 8*boxlength;
vols = zeros(1,100);
rads = (0:100)*bounds/100+0.7;
for kk = 1:100
    vols(1,kk) = 4/3*pi*(rads(1,kk+1)^3-rads(1,kk)^3);
end
raddist_plot= zeros(1,100)

for ll = 1:100
    raddist_plot(1,ll) = 3/2*plotbox(1,ll)/((4*8^3)*vols(1,ll))
end

plot(rads(1,1:100),raddist_plot)
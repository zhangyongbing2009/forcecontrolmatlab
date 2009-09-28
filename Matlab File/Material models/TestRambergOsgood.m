clear all; %close all; clc;
Time = 0.01:0.01:10;
P = 50*sin(Time);

% paprameter for RambergOsgood
% strain = stress/K + A*(stress/K)^(R-1) (original paper)
K = 2.8;  % elstic stiffness
A = 2; %0.5;  % hardening constant
R = 10;   % hardening constant
strain = RambergOsgood_fun(P,K,A,R);

figure;
plot(strain,P)
grid
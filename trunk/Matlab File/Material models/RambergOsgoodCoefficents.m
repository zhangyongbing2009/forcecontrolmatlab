% Calculate the coefficient for the Ramberg-Osgood function
%
% written by T.Y. Yang on 2/27/2009

% clean start
clear all; close all; clc;

% Load the coefficents for the Ramberg-Osgood functions
stress1 = 1.5;
E = 2.8;
strain1 = stress1/E;
stress2 = 1.7; %1.434;
strain2 = 1.96;

% calculate the coefficents for the Ramberg-Osgood functions
minC = [-1000 -1000];
maxC = [1000 1000];
[C,Z] = Coeff_Minimization0(stress1,strain1,E,stress2,strain2,minC,maxC);
A = C(1);
R = C(2);
stress = 0:0.01:stress2;
strain = stress./E.*(1 + A*(stress1./E).^(R-1));
figure;
plot(strain,stress)
hold on
plot(strain1,stress1,'ro')
plot(strain2,stress2,'ro')
grid

% Check material
%
% written by T.Y. Yang on 09/15/2009

% clean start
clear all; close all; clc;

% forcing function
Time = 0.01:0.01:10;
V = 5*sin(Time);
%V = [0.01:0.01:2,2.1,2.1, 1.5:0.01:2.2];

% material property
%Element{1} = 'Hardening';
%Element{1} = 'BiLinearHysteric';
Element{1} = 'NLElastic';
%Element{1} = 'Experimental';
MatData(1).tag    = 1;
MatData(1).E      = 2.8;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.5;    % hardening ratio
MatData(1).Hkin   = 0.1;
MatData(1).amp    = 1;

% initialize the material
feval(Element{1},'initialize',MatData(1));

% loop through the force vector
P = zeros(length(V),1);
for nn = 1:length(P)
    feval(Element{1},'setTrialStrain',MatData(1),V(nn));
    P(nn) = feval(Element{1},'getStress',MatData(1));
    feval(Element{1},'commitState',MatData(1));
end
 
figure;
plot(V,P)
xlabel('Displacement')
ylabel('Force')
grid
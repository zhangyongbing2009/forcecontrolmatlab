% Check material
%
% written by T.Y. Yang on 09/15/2009

% clean start
clear all; close all; clc;

% forcing function
Time = 0.01:0.01:10;
P = 3*sin(Time);
%P = [0.01:0.1:2, 1.5:0.01:2.2];

% material property
%Element{1} = 'HardeningForce';
%Element{1} = 'BiLinearHystericForce';
%Element{1} = 'ExperimentalForce';
Element{1} = 'NLElasticForce';
MatData(1).tag    = 1;
MatData(1).E      = 2.8;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.5;    % hardening ratio
MatData(1).Hkin   = 0.1;
MatData(1).amp   = 5;


% initialize the material
feval(Element{1},'initialize',MatData(1));

% loop through the force vector
V = zeros(length(P),1);
for nn = 1:length(P)
    feval(Element{1},'setTrialStress',MatData(1),P(nn));
    V(nn) = feval(Element{1},'getStrain',MatData(1));
    feval(Element{1},'commitState',MatData(1));
end
 
figure;
plot(V,P)
xlabel('Displacement')
ylabel('Force')
grid
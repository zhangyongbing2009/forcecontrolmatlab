% Nonlinear dynamic analysis using force method
%
% Written by T.Y. Yang and Andreas Schellenberg 09/08/2009

% clean start
clear all;
close all;
clc;

% add the subroutine path to the folder
addpath([pwd '\Material models']);

% element 1 properties - left column
%Element{1} = 'ElasticForce';
%Element{1} = 'BiLinearElasticForce';
%Element{1} = 'HardeningForce';
Element{1} = 'ExperimentalForce';
MatData(1).tag    = 1;
MatData(1).E      = 2.8;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.1;    % hardening ratio
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;

% applied force
T = 1;
Omega = 2*pi/T;
time = 0.01:0.01:T;
Q = .75*sin(time*Omega);
V = zeros(size(Q));
Qm = zeros(size(Q));


% initial flexibility matrix
feval(Element{1},'initialize',MatData(1));

for nn = 1:length(Q);
feval(Element{1},'setTrialStress',MatData(1),Q(nn));
[V(nn),Qm(nn)] = feval(Element{1},'getStrain',MatData(1));
end

% disconnect the connection to xpc
feval(Element{1},'disconnect',MatData(1));
 
% plot the result
figure;
plot(V,Q)
xlabel('V');
ylabel('Q');
grid

figure;
plot(time,Q);
hold on
plot(time,Qm,'r');
grid

figure;
plot(Q,Qm);
grid

% remove the subroutine path to the folder
rmpath([pwd '\Material models']);

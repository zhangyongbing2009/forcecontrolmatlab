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
%Element{1} = 'ExperimentalForce';
Element{1} = 'Experimental';
MatData(1).tag    = 1;
MatData(1).E      = 3.2;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.1;    % hardening ratio
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;

% applied force
T = 1;
Omega = 2*pi/T;
time = 0.01:0.01:T;
V = 0.5*sin(time*Omega);
%V = [0:0.01:2, 1:-0.01:0.0];
Q = zeros(size(V));
Vm = zeros(size(V));

% initial flexibility matrix
feval(Element{1},'initialize',MatData(1));

for nn = 1:length(Q)
   feval(Element{1},'setTrialStrain',MatData(1),V(nn));
   [Q(nn),Vm(nn)] = feval(Element{1},'getStress',MatData(1));
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
plot(Vm,Q,'r')
xlabel('Vm');
ylabel('Q');
grid

figure;
plot(time,V)
hold on
plot(time, Vm,'r')
xlabel('Time');
ylabel('V');
legend('V','Vm')
grid

% remove the subroutine path to the folder
rmpath([pwd '\Material models']);
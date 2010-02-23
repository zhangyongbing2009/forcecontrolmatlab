% Nonlinear dynamic analysis using force method
%
% Written by T.Y. Yang and Andreas Schellenberg 09/08/2009

% clean start
clear all;  close all;  clc;

% add the subroutine path to the folder
addpath([pwd '\Material models']);
addpath([pwd '\Switch schemes']);
addpath([pwd '\Control schemes']);

% element 1 properties
% Element{1} = 'Elastic';
%Element{1} = 'BiLinearElastic';
% Element{1} = 'BiLinearHysteric';
Element{1} = 'Experimental';
MatData(1).tag    = 1;
MatData(1).E      = 2.60;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.02;    % hardening ratio
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;

% applied quantity
T = 1;
Omega = 2*pi/T;
time = 0.01:0.01:3*T;
Damp = 0.15; %in
Famp = Damp*MatData(1).E;
Q = Famp*sin(time*Omega);
V = Damp*sin(time*Omega);
Qm = zeros(size(Q));
Vm = zeros(size(Q));
QIncr = zeros(size(Q));
VIncr = zeros(size(Q));
QmIncr = zeros(size(Q));
VmIncr = zeros(size(Q));

% initial flexibility matrix
feval(Element{1},'initialize',MatData(1));

% % Ramp and Hold and Switch
% feval(Element{1},'setIncrTrialStress',MatData(1),0.5);
% [Vm(1)] = feval(Element{1},'getStrain',MatData(1));
% [Qm(1)] = feval(Element{1},'getStress',MatData(1));
% feval(Element{1},'setIncrTrialStress',MatData(1),0.5);
% [Vm(2)] = feval(Element{1},'getStrain',MatData(1));
% [Qm(2)] = feval(Element{1},'getStress',MatData(1));
% feval(Element{1},'setIncrTrialStress',MatData(1),0.5);
% [Vm(3)] = feval(Element{1},'getStrain',MatData(1));
% [Qm(3)] = feval(Element{1},'getStress',MatData(1));
% %Switch to displacement
% feval(Element{1},'setIncrTrialStrain',MatData(1),0);
% [Vm(4)] = feval(Element{1},'getStrain',MatData(1));
% [Qm(4)] = feval(Element{1},'getStress',MatData(1));
% feval(Element{1},'setIncrTrialStrain',MatData(1),0);
% [Vm(5)] = feval(Element{1},'getStrain',MatData(1));
% [Qm(5)] = feval(Element{1},'getStress',MatData(1));
% feval(Element{1},'setIncrTrialStrain',MatData(1),0);
% [Vm(5)] = feval(Element{1},'getStrain',MatData(1));
% [Qm(5)] = feval(Element{1},'getStress',MatData(1));

for nn = 1:length(Q)/4;
%    QIncr(nn) = Q(nn);
%    feval(Element{1},'setIncrTrialStress',MatData(1),QIncr(nn));
%    [Vm(nn)] = feval(Element{1},'getStrain',MatData(1));
%    [Qm(nn)] = feval(Element{1},'getStress',MatData(1));
   VIncr(nn) = V(nn);
   feval(Element{1},'setIncrTrialStrain',MatData(1),VIncr(nn));
   [Vm(nn)] = feval(Element{1},'getStrain',MatData(1));
   [Qm(nn)] = feval(Element{1},'getStress',MatData(1));
end
feval(Element{1},'commitState',MatData(1));
for nn = length(Q)/4+1:length(Q)/2;   
%    VIncr(nn) = V(nn)-V(length(Q)/4);
%    feval(Element{1},'setIncrTrialStrain',MatData(1),VIncr(nn));
%    [QmIncr(nn)] = feval(Element{1},'getStress',MatData(1));
%    [VmIncr(nn)] = feval(Element{1},'getStrain',MatData(1));
%    Qm(nn) = QmIncr(nn);
%    Vm(nn) = VmIncr(nn);
   QIncr(nn) = Q(nn)-Q(length(Q)/4);
   feval(Element{1},'setIncrTrialStress',MatData(1),QIncr(nn));
   [QmIncr(nn)] = feval(Element{1},'getStress',MatData(1));
   [VmIncr(nn)] = feval(Element{1},'getStrain',MatData(1));
   Qm(nn) = QmIncr(nn);
   Vm(nn) = VmIncr(nn);   
end
feval(Element{1},'commitState',MatData(1));
for nn = length(Q)/2+1:3*length(Q)/4;
%    QIncr(nn) = Q(nn)-Q(length(Q)/2);
%    feval(Element{1},'setIncrTrialStress',MatData(1),QIncr(nn));
%    [VmIncr(nn)] = feval(Element{1},'getStrain',MatData(1));
%    [QmIncr(nn)] = feval(Element{1},'getStress',MatData(1));
%    Vm(nn) = VmIncr(nn);
%    Qm(nn) = QmIncr(nn);
   VIncr(nn) = V(nn)-V(length(Q)/2);
   feval(Element{1},'setIncrTrialStrain',MatData(1),VIncr(nn));
   [VmIncr(nn)] = feval(Element{1},'getStrain',MatData(1));
   [QmIncr(nn)] = feval(Element{1},'getStress',MatData(1));
   Vm(nn) = VmIncr(nn);
   Qm(nn) = QmIncr(nn);
end
feval(Element{1},'commitState',MatData(1));
for nn = 3*length(Q)/4+1:length(Q);
%    VIncr(nn) = V(nn)-V(3*length(Q)/4);
%    feval(Element{1},'setIncrTrialStrain',MatData(1),VIncr(nn));
%    [QmIncr(nn)] = feval(Element{1},'getStress',MatData(1));
%    [VmIncr(nn)] = feval(Element{1},'getStrain',MatData(1));
%    Qm(nn) = QmIncr(nn);
%    Vm(nn) = VmIncr(nn);
   QIncr(nn) = Q(nn)-Q(3*length(Q)/4);
   feval(Element{1},'setIncrTrialStress',MatData(1),QIncr(nn));
   [QmIncr(nn)] = feval(Element{1},'getStress',MatData(1));
   [VmIncr(nn)] = feval(Element{1},'getStrain',MatData(1));
   Qm(nn) = QmIncr(nn);
   Vm(nn) = VmIncr(nn);
end
feval(Element{1},'commitState',MatData(1));
feval(Element{1},'disconnect',MatData(1));
 
% plot the result
figure;
plot(Vm,Qm)
xlabel('Vm');
ylabel('Qm');
grid

figure;
plot(time,Q);
hold on
plot(time,Qm,'r');
plot(time,QIncr,'k');
title('Force Time History')
grid

figure;
plot(time,V);
hold on
plot(time,Vm,'r');
plot(time,VIncr,'k');
title('Displacement Time History')
grid

% figure;
% plot(Q,Qm);
% grid

% remove the subroutine path to the folder
rmpath([pwd '\Control schemes']);
rmpath([pwd '\Switch schemes']);
rmpath([pwd '\Material models']);
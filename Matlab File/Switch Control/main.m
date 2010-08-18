% The main gateway program that runs switch control schemes
%
% Written by Hong Kim
% Created: 10/19/2009
% Last Update: 10/19/09

% clean start
clear all;  close all;  clc;

% add the subroutine path to the folder
% addpath([pwd '\Material models']);
% addpath([pwd '\Switch schemes']);
% addpath([pwd '\Control schemes']);
addpath([pwd '/Material models']);
addpath([pwd '/Switch schemes']);
addpath([pwd '/Control schemes']);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Model
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 Element Model
MODEL = createModel();
% MODEL = createModelT();
% % 4 Element Model
% MODEL = createModel4();

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load GroundMotion Data
%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the ground motion
% GMDir = 'D:\Switch Control\Ground motions\';
GMDir = '/Users/hongkim/Research/Force Control/forcecontrolmatlab/Ground motions/';
dt = 0.02;
SF = 0.15;
g = 386.1;
ag0 = load(fullfile(GMDir,'elcentro.txt'));
t0 = 0:length(ag0)-1;
t0 = dt*t0;
tEnd = t0(end);
ag0 = SF*g*ag0;
% f_RspSpc(ag0,g,dt,0.0,0.01,0.005,3);

% change to analysis deltaT
deltaT = 0.02;
t = deltaT*(0:floor(tEnd/deltaT))';
ag = interp1(t0,ag0,t);
b = [1; 1];
npts = 350;% length(ag);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%
ANALYSIS = initializeAnalysis(deltaT);

% initialize global response variables
U = zeros(MODEL.ndf,npts);
Udot = zeros(MODEL.ndf,npts);
Udotdot = zeros(MODEL.ndf,npts);
Pr = zeros(MODEL.ndf,npts);
errorNorms = zeros(1,npts);
iters = zeros(1,npts);
controlModes = zeros(1,npts);
switchFlags = zeros(1,npts);
ks = zeros(1,npts);
normKs = zeros(1,npts);
normdds = zeros(1,npts);

% initialize element variables
u = zeros(MODEL.numElem,npts);
pr = zeros(MODEL.numElem,npts);
uall = zeros(MODEL.numElem,npts);
prall = zeros(MODEL.numElem,npts);

STATE.prPrev = zeros(MODEL.numElem,1);
STATE.uPrev = zeros(MODEL.numElem,1);
STATE.offsetu = zeros(MODEL.numElem,1);
STATE.offsetpr = zeros(MODEL.numElem,1);
K = MODEL.K;
f = MODEL.f;

% calculations for each time step, i
for i=1:npts-1
   clc
   i
   % get new response quantities
   STATE.U = U(:,i);
   STATE.Udot = Udot(:,i);
   STATE.Udotdot = Udotdot(:,i);
   STATE.Pr = Pr(:,i);
   STATE.pr = pr(:,i);
   STATE.u = u(:,i);
   STATE.i = i;
   STATE.controlMode = controlModes(:,i);
   STATE.switchFlag = switchFlags(:,i);
   MODEL.K = K;
   MODEL.f = f;
   
   % get applied loads
   STATE.Ptp1 = -MODEL.M*b*ag(i+1);
   
   % Run Hybrid Simulation using switch scheme
   [model analysis state] = feval(ANALYSIS.schemeSwitch, MODEL, ANALYSIS, STATE);
   
   STATE.prPrev = pr(:,i);
   STATE.uPrev = u(:,i);
   K = model.K;
   f = model.f;
   U(:,i+1) = state.U;
   Udot(:,i+1) = state.Udot;
   Udotdot(:,i+1) = state.Udotdot;
   Pr(:,i+1) = state.Pr;
   pr(:,i+1) = state.pr;
   u(:,i+1) = state.u;
   controlModes(:,i+1) = state.controlMode;
%    switchFlags(:,i) = state.switchFlag;
%    STATE.offsetu = state.offsetu;
%    STATE.offsetpr = state.offsetpr;
%    ks(:,i+1) = state.k;
   iters(:,i+1) = state.iter;
   errorNorms(:,i) = state.errorNorm;
%    normKs(:,i+1) = state.normK;
%   normdds(:,i+1) = state.normdd;
   
end

% disconnect from experimental sites
for j=1:MODEL.numElem
   if isequal(MODEL.Element{j},'Experimental')
      feval(MODEL.Element{j},'disconnect',MODEL.MatData(j));
   end
end
fclose('all');


%%%%%%%%%%%%%%%%%%%
% Post-Processing
%%%%%%%%%%%%%%%%%%%
% plot the figures
t = t(1:npts);

% % plot the element hysteresis
figure;
for j=1:MODEL.numElem
    subplot(MODEL.numElem,1,j);
    plot(u(j,:),pr(j,:),ANALYSIS.plotFlag);
    xlabel(['u' num2str(j)]);
    ylabel(['pr' num2str(j)]);
    grid
end

% plot the element hysteresis
% figure;
% for j=1:numElem
%     subplot(numElem,1,j);
%     plot(uall(j,1:count-1),prall(j,:),ANALYSIS.plotFlag);
%     xlabel(['uall' num2str(j)]);
%     ylabel(['prall' num2str(j)]);
%     grid
% end

% plot the element force history
figure;
for j=1:MODEL.numElem
    subplot(MODEL.numElem,1,j);
    plot(t,pr(j,:),ANALYSIS.plotFlag);
    xlabel('Time [sec]')
    ylabel(['pr' num2str(j)]);
    grid
end

% plot the element displacement history
figure;
for j=1:MODEL.numElem
    subplot(MODEL.numElem,1,j);
    plot(t,u(j,:),ANALYSIS.plotFlag);
    xlabel('Time [sec]')
    ylabel(['u' num2str(j)]);
    grid
end

% plot the Node response history
figure;
for j=1:MODEL.ndf
    subplot(3,MODEL.ndf,j);
    plot(t,U(j,:),ANALYSIS.plotFlag);
    xlabel('Time [sec]')
    ylabel(['U' num2str(j)]);
    grid    
    subplot(3,MODEL.ndf,j+MODEL.ndf);
    plot(t,Udot(j,:),ANALYSIS.plotFlag);
    xlabel('Time [sec]')
    ylabel(['Udot' num2str(j)]);
    grid
    subplot(3,MODEL.ndf,j+2*MODEL.ndf);
    plot(t,Udotdot(j,:),ANALYSIS.plotFlag);
    xlabel('Time [sec]')
    ylabel(['Udotdot' num2str(j)]);
    grid
end

% figure;
% for j=1:ndf
%     subplot(ndf,1,j);
%     plot(R(j,:));
%     xlabel(['[-]' num2str(j)]);
%     ylabel(['R' num2str(j)]);
%     grid
% end

% error
% figure;
% plot(t,errorNorms)
% ylabel('Error')
% xlabel('Time [sec]')
% grid

% Iterations
figure;
plot(t,iters,ANALYSIS.plotFlag)
ylabel('Iterations/Time Step')
xlabel('Time [sec]')
grid on

% Iterations
% figure;
% plot(normKs)
% hold on
% plot(normdds, 'r')
% plot(1000*(normKs./normdds),'g')
% ylabel('part of jacs')
% xlabel('Time [sec]')
% grid

% Control Mode and secant stiffness
figure;
plot(t,controlModes)
hold on
plot(t,switchFlags,'.k')
grid on
plot(t,ks,'r')
% % plot(10e2*errorNorms,'k')
% % plot(t,100*(normKs./normdds),'g')
% ylabel('[-]')
% xlabel('Time [sec]')
% grid


% experimental element trial displacement
figure;
eD1 = load('ElementDisp1.txt');
plot(eD1,ANALYSIS.plotFlag)
ylabel('trialDisp')
xlabel('Step [-]')
grid

% % experimental element trial force
figure;
eF1 = load('ElementForce1.txt');
plot(eF1,ANALYSIS.plotFlag)
ylabel('trialForce')
xlabel('Step [-]')
grid

% experimental element hysteresis - disp control
figure;
plot(eD1,eF1,ANALYSIS.plotFlag)
ylabel('trialForce')
xlabel('trialDisp')
grid

% % experimental element hysteresis - disp control
% figure;
% plot(eD1(1:end-1),eF1(2:end),ANALYSIS.plotFlag)
% ylabel('trialForce')
% xlabel('trialDisp')
% grid
% 
% % experimental element hysteresis - force control
% figure;
% plot(eD1(2:end),eF1(1:end-1),ANALYSIS.plotFlag)
% ylabel('trialForce')
% xlabel('trialDisp')
% grid

% 
% % plot all experimental force
% figure;
% plot(prall(1,2:end),ANALYSIS.plotFlag)
% ylabel('prall')
% xlabel('Step [-]')
% grid


% ground motion
figure;
plot(t,ag(1:length(t)),ANALYSIS.plotFlag)
ylabel('Ag')
xlabel('Time [sec]')
grid

% remove the subroutine path to the folder
rmpath([pwd '/Control schemes']);
rmpath([pwd '/Switch schemes']);
rmpath([pwd '/Material models']);

% rmpath([pwd '\Control schemes']);
% rmpath([pwd '\Switch schemes']);
% rmpath([pwd '\Material models']);

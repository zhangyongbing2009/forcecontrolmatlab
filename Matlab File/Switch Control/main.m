% The main gateway program that runs switch control schemes
%
% Written by Hong Kim
% Created: 10/19/2009
% Last Update: 10/19/09

% clean start
clear all; % close all;  clc;

% add the subroutine path to the folder
%addpath([pwd '\Material models']);
addpath([pwd '/Material models']);
addpath([pwd '/Switch schemes']);
addpath([pwd '/Control schemes']);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Model
%%%%%%%%%%%%%%%%%%%%%%%%%%
MODEL = createModel();

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load GroundMotion Data
%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the ground motion
GMDir = '/Users/hongkim/Research/Force Control/forcecontrolmatlab/Ground motions/';
dt = 0.02;
SF = 1;
g = 386.1;
ag0 = load(fullfile(GMDir,'elcentro.txt'));
t0 = 0:length(ag0)-1;
t0 = dt*t0;
tEnd = t0(end);
ag0 = SF*g*ag0;

% change to analysis deltaT
deltaT = 0.02;
t = deltaT*(0:floor(tEnd/deltaT))';
ag = interp1(t0,ag0,t);
b = [1; 1];
npts = 400;% length(ag);

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
count = 1;
% initialize element variables
u = zeros(MODEL.numElem,npts);
pr = zeros(MODEL.numElem,npts);
uall = zeros(MODEL.numElem,npts);
prall = zeros(MODEL.numElem,npts);

% calculations for each time step, i
for i=1:npts-1
   
   % get new response quantities
   STATE.U = U(:,i);
   STATE.Udot = Udot(:,i);
   STATE.Udotdot = Udotdot(:,i);
   STATE.Pr = Pr(:,i);
   STATE.pr = pr(:,i);
   STATE.u = u(:,i);
   STATE.i = i;
   
   % get applied loads
   STATE.Ptp1 = -MODEL.M*b*ag(i+1);
   
   % Newton-Raphson algorithm
   [model analysis state] = feval(ANALYSIS.scheme, MODEL, ANALYSIS, STATE);
   
   U(:,i+1) = state.U;
   Udot(:,i+1) = state.Udot;
   Udotdot(:,i+1) = state.Udotdot;
   Pr(:,i+1) = state.Pr;
   pr(:,i+1) = state.pr;
   u(:,i+1) = state.u;
   
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
    plot(t,Udot(j,:),ANALYSIS.plotFlag);
    xlabel('Time [sec]')
    ylabel(['Udot' num2str(j)]);
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

% % error
% figure;
% plot(t,errorNorms)
% ylabel('Error')
% xlabel('Time [sec]')
% grid

% % Iterations
% figure;
% plot(t,iters)
% ylabel('Error')
% xlabel('Time [sec]')
% grid
% 
% % experimental element trial displacement
figure;
eD1 = load('ElementDisp1.txt');
plot(eD1,ANALYSIS.plotFlag)
ylabel('trialDisp')
xlabel('Step [-]')
grid

% % experimental element trial displacement
figure;
eF1 = load('ElementForce1.txt');
plot(eF1,ANALYSIS.plotFlag)
ylabel('trialForce')
xlabel('Step [-]')
grid

% experimental element trial displacement
% figure;
% plot(eD1(1:end-1),eF1(2:end),ANALYSIS.plotFlag)
% ylabel('trialForce')
% xlabel('trialDisp')
% grid

% experimental element trial displacement
figure;
plot(eD1(2:end),eF1(1:end-1),ANALYSIS.plotFlag)
ylabel('trialForce')
xlabel('trialDisp')
grid

% 
% % plot all experimental force
% figure;
% plot(prall(1,2:end),ANALYSIS.plotFlag)
% ylabel('prall')
% xlabel('Step [-]')
% grid


% % ground motion
% figure;
% plot(t,ag(1:length(t)),ANALYSIS.plotFlag)
% ylabel('Ag')
% xlabel('Time [sec]')
% grid

% remove the subroutine path to the folder
rmpath([pwd '/Control schemes']);
rmpath([pwd '/Switch schemes']);
rmpath([pwd '/Material models']);

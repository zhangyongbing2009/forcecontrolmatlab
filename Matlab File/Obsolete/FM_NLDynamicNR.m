% Nonlinear dynamic analysis using force method
%
% Written by T.Y. Yang and Andreas Schellenberg 09/08/2009

% clean start
clear all; close all; clc;

% add the subroutine path to the folder
addpath([pwd '\Material models']);

% constants for the three span truss problem
M = [0.04 0; 0 0.02];
C = 1.01*M;
% B = [1 -1 0; 0 1 -1];
% Bi = [1 1; 0 1; 0 0];
% Bx = [1; 1; 1];

% B = [1 -1 0 -1; 0 1 -1 1];
% Bi = [1 1; 0 1; 0 0; 0 0];
% Bx = [1 0; 1 -1; 1 0; 0 1];

B = [1 -1 0 1; 0 1 -1 0];
Bi = [1 1; 0 1; 0 0; 0 0];
Bx = [1 -1; 1 0; 1 0; 0 1];

% number of global DOF in the model
ndf = size(M,1);
% degree of static indeterminancy
nos = size(Bx,2);

% element 1 properties - left column
Element{1} = 'ElasticForce';
%Element{1} = 'BiLinearElasticForce';
%Element{1} = 'HardeningForce';
%Element{1} = 'ExperimentalForce';
MatData(1).tag    = 1;
MatData(1).E      = 2.66;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.5;    % hardening ratio
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;

% element 2 properties - right column
Element{2} = 'ElasticForce';
MatData(2).tag = 2;
MatData(2).E   = 2.0;

% element 3 properties - spring
Element{3} = 'ElasticForce';
MatData(3).tag = 3;
MatData(3).E   = 5.6;

% element 4 properties - spring
Element{4} = 'ElasticForce';
MatData(4).tag = 4;
MatData(4).E   = 10;

% initial flexibility matrix
numElem = length(Element);
Ft = zeros(numElem);
for j=1:numElem
   feval(Element{j},'initialize',MatData(j));
   Ft(j,j) = feval(Element{j},'getInitialTangent',MatData(j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load GroundMotion Data
%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the ground motion
GMDir = 'D:\Force Control\Ground motions\';
dt = 0.02;
SF = 1;
g = 386.1;
ag0 = load([GMDir 'elcentro.txt']);
t0 = 0:length(ag0)-1;
t0 = dt*t0;
tEnd = t0(end);
ag0 = SF*g*ag0;

% change to analysis deltaT
deltaT = 0.02;
t = deltaT*(0:floor(tEnd/deltaT))';
ag = interp1(t0,ag0,t);
b = [1; 1];
npts = 400; %length(ag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newmark Transient Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis parameters
beta  = 1/4;
gamma = 1/2;
c1 = 1.0;
c2 = gamma/(beta*deltaT);
c3 = 1.0/(beta*deltaT*deltaT);
a1 = (1.0 - gamma/beta);
a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
a3 = -1.0/(beta*deltaT);
a4 = 1.0 - 0.5/beta;

tol = 1.0E-12;
maxIter = 20;

% initialize global response variables
Q = zeros(numElem,npts);
V = zeros(numElem,npts);
U = zeros(ndf,npts);
Udot = zeros(ndf,npts);
Udotdot = zeros(ndf,npts);

% assemble the augmented matrices
Sb = [B; Bx'*Ft];
Mb = [M; zeros(nos,ndf)];
Cb = [C; zeros(nos,ndf)];

% calculations for each time step, nn
for nn = 2:npts
   iter = 1;
   errorNorm = 1;
   Pr = B*Q(:,nn);
   Prb = [Pr; zeros(nos,1)];
   Pb = -Mb*b*ag(nn);
   
   % calculate the initial trial responses
   Q(:,nn) = Q(:,nn-1);
   U(:,nn) = U(:,nn-1);
   Udot(:,nn) = a1*Udot(:,nn-1) + a2*Udotdot(:,nn-1);
   Udotdot(:,nn) = a3*Udot(:,nn-1) + a4*Udotdot(:,nn-1);
   
   while (errorNorm >= tol && iter <= maxIter)
      % set trial forces in elements
      for j=1:numElem
         feval(Element{j},'setTrialStress',MatData(j),Q(j,nn));
      end
      
      % get displacements and flexibilities from elements
      for j=1:numElem
         V(j,nn) = feval(Element{j},'getStrain',MatData(j));
         Ft(j,j) = feval(Element{j},'getTangent',MatData(j));
      end
      
      U(:,nn) = Bi'*V(:,nn);
      Udot(:,nn) = c2*(U(:,nn)-U(:,nn-1)) + a1*Udot(:,nn-1) + a2*Udotdot(:,nn-1);
      Udotdot(:,nn) = c3*(U(:,nn)-U(:,nn-1)) + a3*Udot(:,nn-1) + a4*Udotdot(:,nn-1);
      
      Pr = B*Q(:,nn);
      Prb = [Pr; zeros(nos,1)];
      Sb = [B; Bx'*Ft];
      
      % assemble the Jacobian matrix
      R = Mb*Udotdot(:,nn) + Cb*Udot(:,nn) + Prb - Pb;
      dRdQ = c3*Mb*Bi'*Ft + c2*Cb*Bi'*Ft + Sb;
      dQ = -dRdQ\R;
      % update variables
      Q(:,nn) = Q(:,nn) + dQ;
      
      % update the tolerance and iteration number
      % get energy increment
      errorNorm = norm(dQ);
      iter = iter+1;
   end
   
   if (iter < maxIter)
      % commit the elements
      for j=1:numElem
         feval(Element{j},'commitState',MatData(j));
      end
   else
      error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(nn),...
         ', errorNorm = ',num2str(errorNorm)]);
   end
end

% disconnect from experimental sites
for j=1:numElem
   if isequal(Element{j},'ExperimentalForce')
      feval(Element{j},'disconnect',MatData(j));
   end
end

t = t(1:npts);

figure;
subplot(3,1,1)
plot(V(1,:),Q(1,:));
xlabel('V1');
ylabel('Q1');
grid
subplot(3,1,2)
plot(V(2,:),Q(2,:));
xlabel('V2');
ylabel('Q2');
grid
subplot(3,1,3)
plot(V(3,:),Q(3,:));
xlabel('V3');
ylabel('Q3');
grid

figure;
subplot(3,1,1)
plot(t,Q(1,:))
hold on
ylabel('Q1')
xlabel('Time [sec]')
grid
subplot(3,1,2)
plot(t,Q(2,:))
ylabel('Q2')
xlabel('Time [sec]')
grid
subplot(3,1,3)
plot(t,Q(3,:))
ylabel('Q3')
xlabel('Time [sec]')
grid

figure;
subplot(3,2,1)
plot(t,U(1,:))
ylabel('U1')
grid
subplot(3,2,2)
plot(t,U(2,:))
ylabel('U2')
grid
subplot(3,2,3)
plot(t,Udot(1,:))
ylabel('U1dot')
xlabel('Time [sec]')
grid
subplot(3,2,4)
plot(t,Udot(2,:))
ylabel('U2dot')
grid
subplot(3,2,5)
plot(t,Udotdot(1,:))
ylabel('U1dotdot')
grid
subplot(3,2,6)
plot(t,Udotdot(2,:))
ylabel('U2dotdot')
xlabel('Time [sec]')
grid

% remove the subroutine path to the folder
rmpath([pwd '\Material models']);
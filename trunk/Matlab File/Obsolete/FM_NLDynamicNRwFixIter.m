% Nonlinear dynamic analysis using force method
% Same as DynamicIntegratedForceMethod2.m, except fixed number of iteration is used
%
% Written by T.Y. Yang and Andreas Schellenberg 09/08/2009

% clean start
clear all; %close all; clc;

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
% Element{1} = 'BiLinearElasticForce';
%Element{1} = 'HardeningForce';
%Element{1} = 'ExperimentalForce';
MatData(1).tag    = 1;
MatData(1).E      = 2.66;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.5;    % hardening ratio
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;
MatData(1).Hkin   = 0.01;

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

% number of element 
numElem = size(B,2);

% initial flexibility matrix
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

% fix number of iterations
maxIter = 20;
Tol = 1e10;

% initialize global response variables
Q = zeros(numElem,npts);
Qm = zeros(1,npts);
errorNorms = zeros(1,npts);
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
   Pr = B*Q(:,nn);
   Prb = [Pr; zeros(nos,1)];
   Pb = -Mb*b*ag(nn);
   
   % calculate the initial trial responses
   Q(:,nn) = Q(:,nn-1);
   U(:,nn) = U(:,nn-1);
   Udot(:,nn) = a1*Udot(:,nn-1) + a2*Udotdot(:,nn-1);
   Udotdot(:,nn) = a3*Udot(:,nn-1) + a4*Udotdot(:,nn-1);
   
   for iter=1:maxIter
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
      %Prb = [Pr; zeros(nos,1)];
      Prb = [Pr; Bx'*V(:,nn)];
      Sb = [B; Bx'*Ft];
      
      % assemble the Jacobian matrix
      R = Mb*Udotdot(:,nn) + Cb*Udot(:,nn) + Prb - Pb;
      dRdQ = c3*Mb*Bi'*Ft + c2*Cb*Bi'*Ft + Sb;
      dQ = -dRdQ\R;
      %x = iter/maxIter;
      %scaleddQ = x*(Q(:,nn) + dQ) - (x-1)*Q(:,nn-1) - Q(:,nn);
      scaleddQ = dQ/(maxIter-iter+1);
      
      % update variables
      Q(:,nn) = Q(:,nn) + scaleddQ;
   end
   
   % set trial forces in elements
   for j=1:numElem
      feval(Element{j},'setTrialStress',MatData(j),Q(j,nn));
   end
   
%    get displacements and flexibilities from elements
   for j=1:numElem
      %[V(j,nn),Qm(j,nn)] = feval(Element{j},'getStrain',MatData(j));
      V(j,nn) = feval(Element{j},'getStrain',MatData(j));
   end
   
   U(:,nn) = Bi'*V(:,nn);
   Udot(:,nn) = c2*(U(:,nn)-U(:,nn-1)) + a1*Udot(:,nn-1) + a2*Udotdot(:,nn-1);
   Udotdot(:,nn) = c3*(U(:,nn)-U(:,nn-1)) + a3*Udot(:,nn-1) + a4*Udotdot(:,nn-1);
   
   % commit the elements
   for j=1:numElem
      feval(Element{j},'commitState',MatData(j));
   end
      
   % update the tolerance
   errorNorms(nn) = norm(dQ);
   %errorNorms(nn) = norm(R);
   
   if errorNorms(nn) > Tol
       error(['errorNorms = ' num2str(errorNorms(nn)) ' > Tol = ' num2str(Tol)]);
   end
end

% disconnect from experimental sites
for j=1:numElem
   if isequal(Element{j},'ExperimentalForce')
      feval(Element{j},'disconnect',MatData(j));
   end
end

%%%%%%%%%%%%%%%%%%
% plot the figures
%%%%%%%%%%%%%%%%%%
t = t(1:nn);

% plot the element hysteresis
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(V(j,:),Q(j,:));
    xlabel(['V' num2str(j)]);
    ylabel(['Q' num2str(j)]);
    grid
end

% plot the element force history
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(t,Q(j,:));
    xlabel('Time [sec]')
    ylabel(['Q' num2str(j)]);
    grid
end

% plot the element displacement history
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(t,V(j,:));
    xlabel('Time [sec]')
    ylabel(['V' num2str(j)]);
    grid
end

% plot the Node response history
figure;
for j=1:ndf
    subplot(3,ndf,j);
    plot(t,U(j,:));
    xlabel('Time [sec]')
    ylabel(['U' num2str(j)]);
    grid    
    subplot(3,ndf,j+ndf);
    plot(t,Udot(j,:));
    xlabel('Time [sec]')
    ylabel(['Udot' num2str(j)]);
    grid
    subplot(3,ndf,j+2*ndf);
    plot(t,Udot(j,:));
    xlabel('Time [sec]')
    ylabel(['Udot' num2str(j)]);
    grid
end

% ground motion
figure;
plot(t,ag0(1:length(t)))
ylabel('Ag')
xlabel('Time [sec]')
grid

% remove the subroutine path to the folder
rmpath([pwd '\Material models']);
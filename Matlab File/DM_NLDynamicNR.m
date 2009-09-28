% Nonlinear dynamic analysis using displacement method
% Similar to the DynamicIntegratedForceMethod3.m, except force method is used
%
% Written by T.Y. Yang and Andreas Schellenberg and Hong Kim 09/14/2009

% clean start
clear all; close all; clc;

% add the subroutine path to the folder
addpath([pwd '\Material models']);

% constants for the three span truss problem
M = [0.04 0; 0 0.02];
C = 1.01*M;

% Equilibrium matrix
B = [1 -1 0 1; 0 1 -1 0];
A = B';

% number of global DOF in the model
ndf = size(M,1);

% element 1 properties
%Element{1} = 'Elastic';
%Element{1} = 'BiLinearElastic';
Element{1} = 'BiLinearHysteric';
%Element{1} = 'Hardening';
%Element{1} = 'Experimental';
MatData(1).tag    = 1;
MatData(1).E      = 2.8;
MatData(1).Fy     = 1.5;    % yield stress
MatData(1).b      = 0.1;    % hardening ratio
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;
MatData(1).Hkin   = 0.01;

% element 2 properties
Element{2} = 'Elastic';
MatData(2).tag = 2;
MatData(2).E   = 2.0;

% element 3 properties
Element{3} = 'Elastic';
MatData(3).tag = 3;
MatData(3).E   = 5.6;

% element 4 properties
Element{4} = 'Elastic';
MatData(4).tag = 4;
MatData(4).E   = 10;

% number of element 
numElem = size(B,2);

% initial stiffness matrix
k = zeros(numElem);
for j=1:numElem
   feval(Element{j},'initialize',MatData(j));
   k(j,j) = feval(Element{j},'getInitialTangent',MatData(j));
end
K = A'*k*A;

% calculate natural frequencies and periods
lambda = eig(K,M);
omega  = sort(sqrt(lambda));
T = 2.0*pi./omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load GroundMotion Data
%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the ground motion
GMDir = 'C:\Documents and Settings\Tony\Desktop\Force control\Example\';
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
theta = 0.80;

% initialize global response variables
U = zeros(ndf,npts);
Udot = zeros(ndf,npts);
Udotdot = zeros(ndf,npts);
Pr = zeros(ndf,npts);
errorNorms = zeros(1,npts);

% initialize element variables
u = zeros(numElem,npts);
pr = zeros(numElem,npts);

% set the constants
c1 = 1.0;
c2 = gamma/(beta*dt);
c3 = 1.0/(beta*dt*dt);
a1 = (1.0 - gamma/beta);
a2 = (dt)*(1.0 - 0.5*gamma/beta);
a3 = -1.0/(beta*dt);
a4 = 1.0 - 0.5/beta;

% max iterations and Tol
maxIter = 20;
Tol = 1e-3;

% calculations for each time step, i
for i=1:npts-1
   
   % get new response quantities
   U(:,i+1) = U(:,i);
   Udot(:,i+1) = a1*Udot(:,i) + a2*Udotdot(:,i);
   Udotdot(:,i+1) = a3*Udot(:,i) + a4*Udotdot(:,i);
   
   % get applied loads
   Ptp1 = -M*b*ag(i+1);
   
   % Newton-Raphson algorithm
   iter = 0;
   errorNorm = 1.0;
   while ((errorNorm >= Tol) && (iter <= maxIter))
      
      % transform displacements from global to element DOF
      u(:,i+1) = A*U(:,i+1);
      
      % set trial response in elements
      for j=1:numElem
         feval(Element{j},'setTrialStrain',MatData(j),u(j,i+1));
      end
      
      % get resisting forces and stiffness from elements
      for j=1:numElem
         pr(j,i+1) = feval(Element{j},'getStress',MatData(j));
         k(j,j)    = feval(Element{j},'getTangent',MatData(j));
      end
      
      % transform forces and stiffness from element to global DOF
      Pr(:,i+1) = A'*pr(:,i+1);
      K = A'*k*A;
      
      % get rhs and jacobian
      F  = M*Udotdot(:,i+1) + C*Udot(:,i+1) + Pr(:,i+1) - Ptp1;
      DF = c3*M + c2*C + c1*K;
      
      % solve for displacement increments
      deltaU = DF\(-F);
      
      U(:,i+1) = U(:,i+1) + c1*deltaU;
      Udot(:,i+1) = Udot(:,i+1) +c2*deltaU;
      Udotdot(:,i+1) = Udotdot(:,i+1) + c3*deltaU;
      
      % get energy increment
      errorNorm = norm(F);
      iter = iter + 1;
   end
   
   % commit the elements
   for j=1:numElem
       feval(Element{j},'commitState',MatData(j));
   end
   
   errorNorms(i) = errorNorm;
   
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
   if isequal(Element{j},'Experimental')
      feval(Element{j},'disconnect',MatData(j));
   end
end

% plot the figures
t = t(1:npts);

% plot the element hysteresis
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(u(j,:),pr(j,:));
    xlabel(['u' num2str(j)]);
    ylabel(['pr' num2str(j)]);
    grid
end

% plot the element force history
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(t,pr(j,:));
    xlabel('Time [sec]')
    ylabel(['pr' num2str(j)]);
    grid
end

% plot the element displacement history
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(t,u(j,:));
    xlabel('Time [sec]')
    ylabel(['u' num2str(j)]);
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

% error
figure;
plot(t,errorNorms)
ylabel('Error')
xlabel('Time [sec]')
grid

% ground motion
figure;
plot(t,ag0(1:length(t)))
ylabel('Ag')
xlabel('Time [sec]')
grid

% remove the subroutine path to the folder
rmpath([pwd '\Material models']);
% ONEBAYFRAME_NEWMARK to perform a hybrid simulation of a one-bay-frame
%
% This program interfaces with OpenFresco. The structure is a two
% degree of freedom system as seen in the diagram below.  It uses
% the Newmark time integration scheme.
%
%                        DoF 1     Element 3     DoF 2
%                          [+]o-->--\/\/\/\----o[+]-->
%                           |                    |
%                 Element 1 |                    |  Element 2
%                           |                    |
%                           |                    |
%                        --[+]--              --[+]--
%                        ///////              ///////

% Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
% Created: 03/07
% Revision: A

%close all;
clear all;
clc;

% add the subroutine path to the folder
addpath([pwd '\Material models\']);

%%%%%%%%%%%%%%%%%%%%
% Model Generation
%%%%%%%%%%%%%%%%%%%%
% number of global DOF in the model
ndf = 2;

% define mass matrix
M = zeros(ndf,ndf);
M(1,1) = 0.04;   % mass for DOF 1 (kips/g)
M(2,2) = 0.02;   % mass for DOF 2 (kips/g)

% element 1 properties - left column
Element{1} = 'Elastic';
%Element{1} = 'Hardening';
%Element{1} = 'BiLinearElastic';
MatData(1).tag    = 1;
MatData(1).E      = 2.8;
MatData(1).Fy     = 1.5;
MatData(1).b      = 0.1;
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;
MatData(1).fy     = 1.5;
MatData(1).Hkin   = 0.25/(1.0-0.25)*2.8;
MatData(1).Kiso   = 0.0;

% element 2 properties - right column
Element{2} = 'Elastic';
MatData(2).tag  = 2;
MatData(2).E    = 5.6;
MatData(2).fy   = 500;
MatData(2).Hkin = 0.0;
MatData(2).Kiso = 0.0;

% element 3 properties - spring
Element{3} = 'Elastic';
MatData(3).tag = 3;
MatData(3).E   = 2.0;

% transformation matrix
A = [ 1  0;
      0  1;
     -1  1];

% initial stiffness matrix
numElem = length(Element);
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

% mass-proportional damping matrix
zeta = 0.05;
alphaM = 1.01;  %2.0*zeta*omega(1);
C = alphaM*M;

% Rayleigh damping matrix
%zeta = [0.05, 0.05];
%coef = [1/omega(1) omega(1); 1/omega(2) omega(2)]\[zeta(1); zeta(2)]*2;
%alphaM = coef(1);
%betaKi = coef(2);
%C = alphaM*M + betaKi*K;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load GroundMotion Data
%%%%%%%%%%%%%%%%%%%%%%%%%%
path = 'C:\Documents and Settings\Tony\Desktop\Force control\Models\ThreeSpanTruss\';
dt = 0.02;
g = 386.1;
[t0 ag0] = ReadWriteTHFile('readTHF',[path 'elcentro.txt']);
t0 = dt*t0;
tEnd = t0(end);
ag0 = g*ag0;

% change to analysis deltaT
deltaT = 0.02;
t = deltaT*(0:floor(tEnd/deltaT))';
ag = interp1(t0,ag0,t);
b = [1; 1];
npts = length(ag);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newmark Transient Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis parameters
beta  = 0.25;
gamma = 0.50;
theta = 0.80;
tol = 1.0E-12;
maxIter = 25;

% initialize global response variables
U = zeros(ndf,npts);
Udot = zeros(ndf,npts);
Udotdot = zeros(ndf,npts);
Pr = zeros(ndf,npts);

% initialize element variables
u = zeros(numElem,npts);
pr = zeros(numElem,npts);

% calculations for each time step, i
for i=1:npts-1
   
   % set the constants
   c1 = 1.0;
   c2 = gamma/(beta*dt);
   c3 = 1.0/(beta*dt*dt);
   a1 = (1.0 - gamma/beta);
   a2 = (dt)*(1.0 - 0.5*gamma/beta);
   a3 = -1.0/(beta*dt);
   a4 = 1.0 - 0.5/beta;

   % get new response quantities
   U(:,i+1) = U(:,i);
   Udot(:,i+1) = a1*Udot(:,i) + a2*Udotdot(:,i);
   Udotdot(:,i+1) = a3*Udot(:,i) + a4*Udotdot(:,i);
   
   % get applied loads
   Ptp1 = -M*b*ag(i+1);

   % Newton-Raphson algorithm
   iter = 0;
   energyIncr = 1.0;
   while ((energyIncr >= tol) && (iter <= maxIter))
      
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
      
      % update response quantities
      U(:,i+1) = U(:,i+1) + theta*c1*deltaU;
      Udot(:,i+1) = Udot(:,i+1) + theta*c2*deltaU;
      Udotdot(:,i+1) = Udotdot(:,i+1) + theta*c3*deltaU;
      
      % get energy increment
      energyIncr = 0.5*norm(deltaU'*F);
      iter = iter + 1;
   end
   
   if (iter < maxIter)
      U(:,i+1) = U(:,i+1) + (1.0-theta)*c1*deltaU;
      % commit the elements
      for j=1:numElem
         feval(Element{j},'commitState',MatData(j));
      end
   else
      error('failed to converge in Newton-Raphson algorithm');
   end
end

% disconnect from experimental sites
for j=1:numElem
   if isequal(Element{j},'Experimental')
      feval(Element{j},'disconnect',MatData(j));
   end
end


%%%%%%%%%%%%%%%%%%%
% Post-Processing
%%%%%%%%%%%%%%%%%%%
% plot displacement vs. time
figure;
%figure(1);
plot(t,U(1,:),'b-');
hold('on');
plot(t,U(2,:),'r-');
grid('on');
ylabel('Displacement [in.]');
xlabel('Time [sec]');
legend('DOF 1','DOF 2',1);

% plot element force vs. element deformation
figure;
%figure(2);
plot(u(1,:),pr(1,:),'b-');
hold('on');
plot(u(2,:),pr(2,:),'r-');
grid('on');
ylabel('Resisting Force [kips]');
xlabel('Displacement [in.]');
legend('Element 1','Element 2',2);

figure;
plot(t,pr(1,:));
grid('on');


% Error = M(1,1)*Udotdot(1,:) + C(1,1)*Udot(1,:) + K(1,1)*U(1,:) + K(1,2)*U(2,:) + M(1,1)*ag;
% 
% figure;
% %figure(3);
% plot(t,Error);

% remove the subroutine path to the folder
rmpath([pwd 'Material models\']);

% Nonlinear dynamic analysis using displacement method
%
% Written by T.Y. Yang and Andreas Schellenberg and Hong Kim 09/14/2009

% clean start
clear all;
close all;
clc;

% add the subroutine path to the folder
%addpath([pwd '\Material models']);
addpath([pwd '/Material models']);


%%%%%%%%%%%%%%%%%%%%
% Model Generation
%%%%%%%%%%%%%%%%%%%%
% global mass matrix
M = [0.04 0; 0 0.02];

% equilibrium matrix
B = [1 -1 0; 0 1 -1];
% B = [1 -1 0 -1; 0 1 -1 1];
% B = [1 -1 0 1; 0 1 -1 0];
A = B';

% number of global DOF in the model
ndf = size(M,1);

% element 1 properties
%Element{1} = 'Elastic';
Element{1} = 'BiLinearElastic';
% Element{1} = 'BiLinearHysteretic';
%Element{1} = 'Hardening';
%Element{1} = 'NLElastic';
%Element{1} = 'Experimental';
MatData(1).tag    = 1;
MatData(1).E      = 2.8;
MatData(1).fy     = 1.5;      % yield stress
MatData(1).b      = 1000;     % hardening ratio
MatData(1).Hkin   = MatData(1).b/(1-MatData(1).b)*MatData(1).E;
MatData(1).Kiso   = 0.0;
MatData(1).ipAddr = '127.0.0.1';
MatData(1).ipPort = 8090;

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

% number of elements
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

% mass-proportional damping matrix
zeta = 0.05;
alphaM = 2.0*zeta*omega(1);
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
% load the ground motion
GMDir = pwd;
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
npts = 100;% length(ag);


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

% max iterations and tol
maxIter = 100;
tol = 1.0E-6;

% initialize global response variables
U = zeros(ndf,npts);
Udot = zeros(ndf,npts);
Udotdot = zeros(ndf,npts);
Pr = zeros(ndf,npts);
errorNorms = zeros(1,npts);
iters = zeros(1,npts);
count = 1;
% initialize element variables
u = zeros(numElem,npts);
pr = zeros(numElem,npts);
uall = zeros(numElem,npts);
prall = zeros(numElem,npts);

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
   while ((errorNorm >= tol) && (iter <= maxIter))
      
      % get resisting forces and stiffness from elements
      for j=1:numElem
         pr(j,i+1) = feval(Element{j},'getStress',MatData(j));
         k(j,j)    = feval(Element{j},'getTangent',MatData(j));
      end
      
      prall(:,count) = pr(:,i+1);
      % transform forces and stiffness from element to global DOF
      Pr(:,i+1) = A'*pr(:,i+1);
      K = A'*k*A;
      
      % get rhs and jacobian
      R(:,count)  = M*Udotdot(:,i+1) + C*Udot(:,i+1) + Pr(:,i+1) - Ptp1;
      dRdU = c3*M + c2*C + c1*K;
      
      % solve for displacement increments
      deltaU = dRdU\(-R(:,count));
      
      % update response quantities
      U(:,i+1) = U(:,i+1) + c1*deltaU;
      Udot(:,i+1) = Udot(:,i+1) +c2*deltaU;
      Udotdot(:,i+1) = Udotdot(:,i+1) + c3*deltaU;
      
      % transform displacements from global to element DOF
      u(:,i+1) = A*U(:,i+1);
      uall(:,count+1) = u(:,i+1);
      % set trial response in elements
      for j=1:numElem
         feval(Element{j},'setTrialStrain',MatData(j),u(j,i+1));
      end
      
      % update the error norm and iteration number
      errorNorm = norm(deltaU);
      %errorNorm = norm(R);
      iter = iter + 1;
      count = count + 1;
   end
   
   iters(i) = iter;
   errorNorms(i) = errorNorm;
   
   if (iter < maxIter)
       % commit the elements
       for j=1:numElem
           feval(Element{j},'commitState',MatData(j));
       end
   else
       error(['failed to converge in Newton-Raphson algorithm: Step = ',num2str(i),...
           ', errorNorm = ',num2str(errorNorm)]);
   end
   
end

% disconnect from experimental sites
for j=1:numElem
   if isequal(Element{j},'Experimental')
      feval(Element{j},'disconnect',MatData(j));
   end
end
fclose('all');


%%%%%%%%%%%%%%%%%%%
% Post-Processing
%%%%%%%%%%%%%%%%%%%
% plot the figures
t = t(1:npts);

% % plot the element hysteresis
% figure;
% for j=1:numElem
%     subplot(numElem,1,j);
%     plot(u(j,:),pr(j,:),'b-');
%     xlabel(['u' num2str(j)]);
%     ylabel(['pr' num2str(j)]);
%     grid
% end

% plot the element hysteresis
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(uall(j,1:count-1),prall(j,:),'b-');
    xlabel(['uall' num2str(j)]);
    ylabel(['prall' num2str(j)]);
    grid
end

% % plot the element force history
% figure;
% for j=1:numElem
%     subplot(numElem,1,j);
%     plot(t,pr(j,:),'b-');
%     xlabel('Time [sec]')
%     ylabel(['pr' num2str(j)]);
%     grid
% end
% 
% % plot the element displacement history
% figure;
% for j=1:numElem
%     subplot(numElem,1,j);
%     plot(t,u(j,:),'b-');
%     xlabel('Time [sec]')
%     ylabel(['u' num2str(j)]);
%     grid
% end
% 
% % plot the Node response history
% figure;
% for j=1:ndf
%     subplot(3,ndf,j);
%     plot(t,U(j,:),'b-');
%     xlabel('Time [sec]')
%     ylabel(['U' num2str(j)]);
%     grid    
%     subplot(3,ndf,j+ndf);
%     plot(t,Udot(j,:),'b-');
%     xlabel('Time [sec]')
%     ylabel(['Udot' num2str(j)]);
%     grid
%     subplot(3,ndf,j+2*ndf);
%     plot(t,Udot(j,:),'b-');
%     xlabel('Time [sec]')
%     ylabel(['Udot' num2str(j)]);
%     grid
% end

figure;
for j=1:ndf
    subplot(ndf,1,j);
    plot(R(j,:));
    xlabel(['[-]' num2str(j)]);
    ylabel(['R' num2str(j)]);
    grid
end

% % error
% figure;
% plot(t,errorNorms)
% ylabel('Error')
% xlabel('Time [sec]')
% grid

% Iterations
figure;
plot(t,iters)
ylabel('Error')
xlabel('Time [sec]')
grid

% experimental element trial displacement
figure;
eD1 = load('ElementDisp1.txt');
plot(eD1,'b-')
ylabel('trialDisp')
xlabel('Step [-]')
grid

% plot all experimental force
figure;
plot(prall(1,2:end),'b-')
ylabel('prall')
xlabel('Step [-]')
grid


% % ground motion
% figure;
% plot(t,ag(1:length(t)),'b-')
% ylabel('Ag')
% xlabel('Time [sec]')
% grid

% remove the subroutine path to the folder
% rmpath([pwd '\Material models']);
rmpath([pwd '/Material models']);

% Nonlinear dynamic analysis using force method with fixed number of iterations
%
% Written by T.Y. Yang and Andreas Schellenberg 09/08/2009

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

% equilibrium matrices
B = [1 -1 0; 0 1 -1];
Bi = [1 1; 0 1; 0 0];
Bx = [1; 1; 1];

% B = [1 -1 0 -1; 0 1 -1 1];
% Bi = [1 1; 0 1; 0 0; 0 0];
% Bx = [1 0; 1 -1; 1 0; 0 1];

% B = [1 -1 0 1; 0 1 -1 0];
% Bi = [1 1; 0 1; 0 0; 0 0];
% Bx = [1 -1; 1 0; 1 0; 0 1];

% number of global DOF in the model
ndf = size(M,1);
% degree of static indeterminancy
nos = size(Bx,2);

% element 1 properties - left column
%Element{1} = 'ElasticForce';
%Element{1} = 'BiLinearElasticForce';
Element{1} = 'BiLinearHystereticForce';
%Element{1} = 'HardeningForce';
%Element{1} = 'ExperimentalForce';
MatData(1).tag    = 1;
MatData(1).E      = 2.8;
MatData(1).fy     = 1.5;      % yield stress
MatData(1).b      = 0.01;     % hardening ratio
MatData(1).Hkin   = MatData(1).b/(1-MatData(1).b)*MatData(1).E;
MatData(1).Kiso   = 0.0;
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

% number of elements
numElem = size(B,2);

% initial flexibility and global stiffness matrix
f = zeros(numElem);
for j=1:numElem
   feval(Element{j},'initialize',MatData(j));
   f(j,j) = feval(Element{j},'getInitialTangent',MatData(j));
end
K = B*(f\B');

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

% max iterations
maxIter = 20;

% initialize global response variables
U = zeros(ndf,npts);
Udot = zeros(ndf,npts);
Udotdot = zeros(ndf,npts);
Pr = zeros(ndf,npts);
errorNorms = zeros(1,npts);

% initialize element variables
u = zeros(numElem,npts);
pr = zeros(numElem,npts);

% assemble the augmented matrices
Mb = [M; zeros(nos,ndf)];
Cb = [C; zeros(nos,ndf)];

% calculations for each time step, i
for i = 1:npts-1
   
   % get new response quantity
   pr(:,i+1) = pr(:,i);
   
   % get applied loads
   Pb = -Mb*b*ag(i+1);
   
   for iter=1:maxIter
      
      % get displacements and flexibilities from elements
      for j=1:numElem
         u(j,i+1) = feval(Element{j},'getStrain',MatData(j));
         f(j,j) = feval(Element{j},'getTangent',MatData(j));
      end
      
      % update response quantities
      U(:,i+1) = Bi'*u(:,i+1);
      Udot(:,i+1) = c2*(U(:,i+1)-U(:,i)) + a1*Udot(:,i) + a2*Udotdot(:,i);
      Udotdot(:,i+1) = c3*(U(:,i+1)-U(:,i)) + a3*Udot(:,i) + a4*Udotdot(:,i);
      
      % transform forces from element to global DOF
      Pr(:,i+1) = B*pr(:,i+1);
      Prb = [Pr(:,i+1); Bx'*u(:,i+1)];
      Sb = [B; Bx'*f];
      
      % get rhs and jacobian
      R = Mb*Udotdot(:,i+1) + Cb*Udot(:,i+1) + Prb - Pb;
      dRdQ = (c3*Mb + c2*Cb)*Bi'*f + Sb;
      
      % solve for force increments
      deltaQ = dRdQ\(-R);
      
      % substeps
      x = iter/maxIter;
      scaleddeltaQ = x*(pr(:,i+1) + deltaQ) - (x-1)*pr(:,i) - pr(:,i+1);
%       scaleddeltaQ = deltaQ/(maxIter-iter+1);
      
      % update variables
      pr(:,i+1) = pr(:,i+1) + scaleddeltaQ;
      
      % set trial forces in elements
      for j=1:numElem
         feval(Element{j},'setTrialStress',MatData(j),pr(j,i+1));
      end
   end
   
   % get displacements and flexibilities from elements
   for j=1:numElem
      u(j,i+1) = feval(Element{j},'getStrain',MatData(j));
      f(j,j) = feval(Element{j},'getTangent',MatData(j));
   end
   
   % update response quantities
   U(:,i+1) = Bi'*u(:,i+1);
   Udot(:,i+1) = c2*(U(:,i+1)-U(:,i)) + a1*Udot(:,i) + a2*Udotdot(:,i);
   Udotdot(:,i+1) = c3*(U(:,i+1)-U(:,i)) + a3*Udot(:,i) + a4*Udotdot(:,i);
   
   errorNorms(i) = norm(scaleddeltaQ);
   
   % commit the elements
   for j=1:numElem
       feval(Element{j},'commitState',MatData(j));
   end
end

% disconnect from experimental sites
for j=1:numElem
   if isequal(Element{j},'ExperimentalForce')
      feval(Element{j},'disconnect',MatData(j));
   end
end
fclose('all');


%%%%%%%%%%%%%%%%%%%
% Post-Processing
%%%%%%%%%%%%%%%%%%%
% plot the figures
t = t(1:npts);

% plot the element hysteresis
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(u(j,:),pr(j,:),'r-');
    xlabel(['u' num2str(j)]);
    ylabel(['pr' num2str(j)]);
    grid
end

% plot the element force history
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(t,pr(j,:),'r-');
    xlabel('Time [sec]')
    ylabel(['pr' num2str(j)]);
    grid
end

% plot the element displacement history
figure;
for j=1:numElem
    subplot(numElem,1,j);
    plot(t,u(j,:),'r-');
    xlabel('Time [sec]')
    ylabel(['u' num2str(j)]);
    grid
end

% plot the Node response history
figure;
for j=1:ndf
    subplot(3,ndf,j);
    plot(t,U(j,:),'r-');
    xlabel('Time [sec]')
    ylabel(['U' num2str(j)]);
    grid    
    subplot(3,ndf,j+ndf);
    plot(t,Udot(j,:),'r-');
    xlabel('Time [sec]')
    ylabel(['Udot' num2str(j)]);
    grid
    subplot(3,ndf,j+2*ndf);
    plot(t,Udot(j,:),'r-');
    xlabel('Time [sec]')
    ylabel(['Udot' num2str(j)]);
    grid
end

% error
figure;
plot(t,errorNorms,'r-')
ylabel('Error')
xlabel('Time [sec]')
grid

% experimental element trial force
figure;
eF1 = load('ElementForce1.txt');
plot(eF1,'r-')
ylabel('trialForce')
xlabel('Step [-]')
grid

% % ground motion
% figure;
% plot(t,ag(1:length(t)),'r-')
% ylabel('Ag')
% xlabel('Time [sec]')
% grid

% remove the subroutine path to the folder
% rmpath([pwd '\Material models']);
rmpath([pwd '/Material models']);
